#include "InputOutput.h"
#include <fstream>

IO::Writer::Writer(const std::string& path){
    output_file = std::ofstream(path);
    if (!output_file.is_open()) {
        failed = true;
        throw std::string("Unable to open output file at: '" + path + "'.");
    }
    output_file << "id\tmass\tx\ty\tvel_x\tvel_y\n";
}

IO::Writer::~Writer(void){
    if (!failed)
        output_file.close();
}

void IO::Writer::write_body(const IO::Body& io_body){
    if (failed)
        throw std::string("Writer failed.");
    output_file << io_body.id << '\t'
        << io_body.mass << '\t'
        << io_body.x << '\t'
        << io_body.y << '\t'
        << io_body.vel_x << '\t'
        << io_body.vel_y << std::endl;
}