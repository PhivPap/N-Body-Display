#ifndef INPUT_OUTPUT
#define INPUT_OUTPUT

#include <fstream>
#include <string>
#include <vector>

namespace IO {
    class Body {
    public:
        uint32_t id;
        double mass, x, y, vel_x, vel_y;
        Body(void);
        Body(uint32_t id, double mass, double x, double y, double vel_x, double vel_y);
        void copy_from(const Body& b);
    };


    class Parser {
    public:
        Parser(void) = delete;
        Parser(const std::string& path);
        ~Parser(void);
        int next_body_info(Body& io_body);

    private:
        bool failed = false;
        uint32_t parsed_bodies = 0;
        std::ifstream input_file;
    };


    class Writer {
    public:
        Writer(const std::string& path);
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
        int32_t get_next_idx(const std::string& str);  // returns - 1 if not found and -2 if next does not exist
        std::string get(uint32_t idx);
    private:
        const int argc;
        const char** argv;
    };
}

#endif