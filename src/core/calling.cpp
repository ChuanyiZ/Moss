//
// Created by Chuanyi Zhang on 2019-05-23.
//

#include "calling.h"
#include <cmath>
#include <map>
#include <cassert>
#include <cstring>
#include <limits>

using namespace moss;

// TODO: limited number of tumor samples?
// - Use custom arbitrary long binary indicator?
SnvCaller::SnvCaller(int n_tumor_sample, std::string normal, double mu, int max_depth, int grid_size)
    : n_tumor_sample(n_tumor_sample),
      normal_result(std::move(VcfReader(normal))),
      mu(mu),
      max_depth(max_depth),
      gridSize(grid_size),
      is_empty(n_tumor_sample),
      n_tumor(n_tumor_sample * 3),
      n_normal(n_tumor_sample),
      p_err(n_tumor_sample, std::vector<double>(max_depth)),
      is_normal(n_tumor_sample, std::vector<char>(max_depth)),
      is_tumor(n_tumor_sample, std::vector<char>(max_depth * 3)),
      likelihoods(n_tumor_sample, std::vector<std::vector<double>>(3, std::vector<double>(grid_size))) {
    stepSize = 1.0 / (gridSize - 1);
    logNoisePriorComplement = log(1 - mu);
    logPriorZComplement = new double[n_tumor_sample];
    for (int idx = 0; idx < n_tumor_sample; idx++) {
        logPriorZComplement[idx] = logNoisePriorComplement - log((1 << (idx + 1)) - 1);
    }
    logMu = log1p(mu - 1);
    logUniform = -log(gridSize - 1);
    eps = 0.1;
}

SnvCaller::~SnvCaller() {
    delete[](logPriorZComplement);
}

BaseSet SnvCaller::normal_calling(const std::vector<Read> &column, uint8_t ref) {
    std::map<uint8_t, int> count;
    for (const auto &r : column) {
        count[r.base]++;
    }
    for (const auto &c : count) {
        if ((c.first & ref) != 0) {
            if (abs(c.second / static_cast<double>(column.size()) - 0.5) <= eps) {
                return BaseSet(c.first | ref);
            } else if (1 - eps <= c.second / static_cast<double>(column.size())) {
                return BaseSet(c.first);
            }
        }
    }
    return BaseSet(ref);
}

BaseSet SnvCaller::normal_calling(const std::string &contig, locus_t pos, uint8_t ref) {
    RecData search = normal_result.find(contig, pos);
    if (search.bases.is_valid() && search.is_pass) {
        return search.bases;
    } else {
        return BaseSet(ref);
    }
}

void
SnvCaller::calc_likelihood(const std::vector<std::vector<Read>> &aligned, BaseSet normal_bases, BaseSet tumor_base) {
    // pre-calculate
    auto n_gt = tumor_base.size();
    std::fill(n_normal.begin(), n_normal.end(), 0);
    std::fill(n_tumor.begin(), n_tumor.end(), 0);
    int idx_sample = 0;
    for (const auto &sample : aligned) {
        int idx_read = 0;
        for (auto r = sample.begin(); r != sample.end(); ++r, ++idx_read) {
            if (idx_read >= max_depth) {
                max_depth *= 2;
                for (int i = 0; i < n_tumor_sample; ++i) {
                    p_err[i].resize(max_depth);
                    is_normal[i].resize(max_depth);
                    is_tumor[i].resize(max_depth * 3);
                }
            }
            p_err[idx_sample].at(idx_read) = qphred2prob(r->qual);
            is_normal[idx_sample].at(idx_read) = normal_bases.contain(r->base);
            n_normal[idx_sample] += is_normal[idx_sample][idx_read] ? 1 : 0;
            int idx_base = 0;
            for (const auto &tumorBase : tumor_base) {
                bool eq = r->base == tumorBase;
                is_tumor[idx_sample].at(idx_read * n_gt + idx_base) = eq;
                n_tumor[idx_sample * n_gt + idx_base] += eq ? 1 : 0;
                idx_base++;
            }
        }
        idx_sample++;
    }
    // n_tumor_sample x n_tumor_base x gridSize
    idx_sample = 0;
    for (const auto &sample : aligned) {
        int idx_base = 0;
        for (const auto &base : tumor_base) {
            for (int idx_step = 0; idx_step < gridSize; ++idx_step) {
                double lhood = 0;
                double f = idx_step * stepSize;
                unsigned long sample_size = sample.size();
                assert(sample_size < max_depth);
                for (int j = 0; j < sample_size; ++j) {
                    if (is_normal[idx_sample][j]) {
                        lhood += log(f * p_err[idx_sample][j] / 3 +
                                     (1 - f) * (1 + p_err[idx_sample][j] * (double(normal_bases.size()) - 4.0) / 3));
                        assert(lhood <= 0);
                    } else if (is_tumor[idx_sample][j * n_gt + idx_base]) {
                        lhood += log(f * (1 - p_err[idx_sample][j]) + (1 - f) * p_err[idx_sample][j] / 3);
                        assert(lhood <= 0);
                    } else {
                        // f * err / 3 + (1-f) * err / 3
                        lhood += log(p_err[idx_sample][j]) - log(3);
                        assert(lhood <= 0);
                    }
                }
                unsigned int n_t = n_tumor[idx_sample * n_gt + idx_base];
                unsigned long n_n = n_normal[idx_sample];
                double coeff = log_trinomial(sample_size - n_n - n_t, n_n, n_t);
                assert(!std::isinf(coeff));
                likelihoods[idx_sample][idx_base][idx_step] = lhood + coeff;
            }
            ++idx_base;
        }
        ++idx_sample;
    }
}

double
SnvCaller::calling(const std::string &chrom, locus_t pos, const Pileups &pile, BaseSet &normal_gt, uint8_t &tumor_gt,
                   unsigned long &Z, Annotation &annos) {
    const std::vector<std::vector<Read>> &columns = pile.get_read_columns();
    uint8_t ref = pile.get_ref();
    if (normal_result.empty()) {
        normal_gt = normal_calling(columns[0], ref);
    } else {
        normal_gt = normal_calling(chrom, pos, ref);
    }
    BaseSet tumor_bases = normal_gt.complement();
    calc_likelihood(std::vector<std::vector<moss::Read>>(columns.begin() + 1, columns.end()), normal_gt, tumor_bases);
    size_t n_valid_tumor_sample = 0;
    for (int i = 1; i < columns.size(); i++) {
        if (columns[i].size() == 0) {
            is_empty[i - 1] = true;
        } else {
            is_empty[i - 1] = false;
            n_valid_tumor_sample++;
        }
        annos.cnt_read[i - 1] = columns[i].size();
    }

    // pre-compute integral of log likelihood under z = 0, 1
    // integral over frequency, P(D_i | Z_i=z_i, Gn)
    std::vector<std::vector<double> > llh_integral(n_tumor_sample);
    int idx_sample = 0;
    for (const auto &sample : likelihoods) {
        llh_integral[idx_sample].resize(tumor_bases.size());
        for (int idx_base = 0; idx_base < tumor_bases.size(); ++idx_base) {
            const auto &sample_base = sample[idx_base];
            llh_integral[idx_sample][idx_base] = log_sum_exp(sample_base) + logUniform;
        }
        idx_sample++;
    }
    // sum over 2^m-1 of z, and tumor nucleotide
    double max_nume_elem,
        max_deno_elem,
        max_tumor_evidence = -std::numeric_limits<double>::infinity(),
        max_evidence_elem;
    double nume,
        deno;
    log_sum_exp_init(max_nume_elem, nume);
    log_sum_exp_init(max_deno_elem, deno);
    int tumor_gt_idx;
    int idx_nuc = 0;
    for (const auto &tumor_base : tumor_bases) {
        double evidence;
        log_sum_exp_init(max_evidence_elem, evidence);
        for (int z = 0; z < (1 << n_tumor_sample); ++z) {
            double llh = 0;
            for (int idx_sample = 0; idx_sample < n_tumor_sample; ++idx_sample) {
                if (is_empty[idx_sample]) {
                    continue;
                }
                int z_sample = (z >> idx_sample) & 1;
                llh += z_sample ? llh_integral[idx_sample][idx_nuc] : likelihoods[idx_sample][idx_nuc][0];
            }
            if (z == 0) {
                double temp = llh + logMu;
                log_sum_exp_iter(max_nume_elem, nume, temp);
                log_sum_exp_iter(max_evidence_elem, evidence, temp);
            } else {
                log_sum_exp_iter(max_evidence_elem, evidence, llh + logPriorZComplement[n_valid_tumor_sample - 1]);
            }
        }
        evidence = log_sum_exp_final(max_evidence_elem, evidence);
        if (max_tumor_evidence <= evidence) {
            max_tumor_evidence = evidence;
            tumor_gt = tumor_base;
            tumor_gt_idx = idx_nuc;
        }
        log_sum_exp_iter(max_deno_elem, deno, evidence);
        idx_nuc++;
    }
    Z = 0;
    for (int i = 0; i < n_tumor_sample; ++i) {
        if (!is_empty[i] && (likelihoods[i][tumor_gt_idx][0] < llh_integral[i][tumor_gt_idx])) {
            Z |= (1 << i);
            annos.genotype[i] = 1;
        } else {
            annos.genotype[i] = 0;
        }
        annos.cnt_tumor[i] = n_tumor[i * tumor_bases.size() + tumor_gt_idx];
    }

    double log_prob_non_soma = log_sum_exp_final(max_nume_elem, nume) - log_sum_exp_final(max_deno_elem, deno);
    return log_prob_non_soma;
}

double moss::qphred2prob(int qphred) {
    return pow(10.f, -static_cast<double>(qphred) / 10.f);
}

double moss::binom(unsigned int n, unsigned int k) {
    assert(n >= k);
    if (n == 0 || n == 1 || k == 0 || k == n) {
        return 1;
    }
    if (k == 1) {
        return n;
    }
    if (2 * k > n) {
        k = n - k;
    }
    double bin = n - k + 1;
    for (int i = 2, j = n - k + 2; i <= k; ++i, ++j) {
        bin *= j;
        bin /= i;
    }
    return bin;
}

double moss::trinomial(unsigned long s, unsigned long k, unsigned long t) {
    return binom(s + k + t, t) * binom(s + k, k);
}

double moss::log_trinomial(unsigned long s, unsigned long k, unsigned long t) {
    unsigned long sum = s + k + t;
    if (sum == 0) {
        return 0;
    }
    double log_tri;
    int non_zero_count = int(s != 0) + int(k != 0) + int(t != 0);
    if (sum > 651) {
        // Stirling's approx.
        // C++ standard double < 10^308
        // multinomial(217, 217, 217) is the largest without overflow, => 3*217 = 651
        log_tri = (sum + 0.5) * log(sum);
        log_tri -= (s == 0) ? 0 : (s + 0.5) * log(s);
        log_tri -= (k == 0) ? 0 : (k + 0.5) * log(k);
        log_tri -= (t == 0) ? 0 : (t + 0.5) * log(t);
        log_tri -= log(2 * M_PI) * (non_zero_count - 1) / 2;
    } else {
        log_tri = log(trinomial(s, k, t));
    }
    return log_tri;
}

template<typename T>
T moss::log_sum_exp(std::vector<T> array) {
    T max_elem = -std::numeric_limits<T>::infinity(),
        accum{};
    for (const auto &item : array) {
        if (item >= max_elem) {
            accum *= exp(max_elem - item);
            accum += 1.f;
            max_elem = item;
        } else {
            accum += exp(item - max_elem);
        }
    }
    return max_elem + log(accum);
}

template<typename T>
void moss::log_sum_exp_init(T &max_elem, T &accum) {
    max_elem = -std::numeric_limits<T>::infinity();
    accum = T{};
}

template<typename T>
void moss::log_sum_exp_iter(T &max_elem, T &accum, T item) {
    if (item >= max_elem) {
        accum *= exp(max_elem - item);
        accum += 1.f;
        max_elem = item;
    } else {
        accum += exp(item - max_elem);
    }
}

template<typename T>
T moss::log_sum_exp_final(T &max_elem, T &accum) {
    return max_elem + log(accum);
}
