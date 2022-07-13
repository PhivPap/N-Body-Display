#include "InputOutput.h"
#include "../MyTypes.h"

IO::ArgParser::ArgParser(int argc, const char** argv) : argc(argc), argv(argv) {}

i32 IO::ArgParser::get_next_idx(const sstr& str){
    i32 idx = 0;
    for (i32 i = 0; i < argc; i++)
        if (str == argv[i])
            if (i + 1 >= argc)
                return -2;
            else 
                return i + 1;
    return -1;
}

sstr IO::ArgParser::get(u32 idx){
    if (idx >= argc)
        throw sstr("ArgParser error: index out of bounds.");
    return sstr(argv[idx]);
}