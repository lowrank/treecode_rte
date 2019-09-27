/*
 * Configuration file input
 *
 * Yimin Zhong @ yzhong@utexas.edu.
 *
 */
#ifndef LEVELSET_CONFIG_H
#define LEVELSET_CONFIG_H

#include "../utils.h"
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>


class config {
public:
    config();
    ~config();
    config(char* path);
    std::map<std::string, std::string> options;
    void parse(std::istream &cfgFile);
    void print();
};


#endif //LEVELSET_CONFIG_H
