#include "InputOutput.h"


IO::ArgParser::ArgParser(int argc, const char** argv) : argc(argc), argv(argv) {}

int32_t IO::ArgParser::get_next_idx(const std::string& str){
    int32_t idx = 0;
    for (uint32_t i = 0; i < argc; i++)
        if (str == argv[i])
            if (i + 1 >= argc)
                return -2;
            else 
                return i + 1;
    return -1;
}

std::string IO::ArgParser::get(uint32_t idx){
    if (idx >= argc)
        throw std::string("ArgParser error: index out of bounds.");
    return std::string(argv[idx]);
}