#ifndef INPUT_OUTPUT
#define INPUT_OUTPUT

#include <fstream>
#include <string>
#include <vector>
#include "../MyTypes.h"

namespace IO {
    class Body {
    public:
        u32 id;
        f64 mass, x, y, vel_x, vel_y;
        Body(void);
        Body(u32 id, f64 mass, f64 x, f64 y, f64 vel_x, f64 vel_y);
        void copy_from(const Body& b);
    };


    class Parser {
    public:
        Parser(void) = delete;
        Parser(const sstr& path);
        ~Parser(void);
        int next_body_info(Body& io_body);

    private:
        bool failed = false;
        u32 parsed_bodies = 0;
        std::ifstream input_file;
    };


    class Writer {
    public:
        Writer(const sstr& path);
        ~Writer(void);
        void write_body(const Body& io_body);

    private:
        bool failed = false;
        std::ofstream output_file;
    };

    class ArgParser {
    public:
        ArgParser(void) = delete;
        ArgParser(int argc, const char** argv);
        i32 get_next_idx(const sstr& str);  // returns - 1 if not found and -2 if next does not exist
        sstr get(u32 idx);
    private:
        const int argc;
        const char** argv;
    };
}

#endif