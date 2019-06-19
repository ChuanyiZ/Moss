//
// Created by Chuanyi Zhang on 2019-06-03.
//

#include <fstream>
#include <sstream>
#include <vector>
#include "loci.h"

moss::MapContigLoci moss::merge_loci(std::vector<std::string> filenames) {
    MapContigLoci loci;
    for (const auto &filename : filenames) {
        std::string target,
                position;
        std::ifstream file(filename);
        std::string line;
        while (std::getline(file, line, '\t')) {
            auto temp = std::stringstream(line);
            std::getline(temp, target, ':');
            std::getline(temp, position, ':');
            auto item = loci.find(target);
            if (item != loci.end()) {
                item->second.insert(std::stoul(position));
            } else {
                loci.insert(std::make_pair(target, std::set<unsigned long>({std::stoul(position)})));
            }
        }
        file.close();
    }
    return loci;
}
