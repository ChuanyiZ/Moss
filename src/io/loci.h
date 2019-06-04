//
// Created by Chuanyi Zhang on 2019-06-03.
//

#ifndef MOSS_LOCI_H
#define MOSS_LOCI_H

#include <utility>
#include <map>
#include <set>
#include <string>

namespace moss {
    std::map<std::string, std::set<unsigned long>> merge_loci(std::vector<std::string> filenames);
}

#endif //MOSS_LOCI_H
