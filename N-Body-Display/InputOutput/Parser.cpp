#include "InputOutput.h"
#include <fstream>
#include <vector>
#include <sstream>
#include <iostream>


IO::Parser::Parser(const std::string& path){
    std::string line, id;
    input_file = std::ifstream(path);

    if (!input_file.is_open()) {
        failed = true;
        throw std::string("Unable to open input file at: '" + path + "'.");
    }

    if (!getline(input_file, line)){
        failed = true;
        throw std::string("Parser received empty file.");
    }
}

IO::Parser::~Parser(void){
    if (!failed)
        input_file.close();
}

// code from: https://stackoverflow.com/questions/1120140/how-can-i-read-and-parse-csv-files-in-c
static std::vector<std::string> tokenize_string(const std::string& str, char delim){
    std::vector<std::string> tokens;
    std::stringstream str_stream(str);
    std::string token;
    while (std::getline(str_stream, token, delim))
        tokens.push_back(token);
    if (!str_stream && token.empty())
        tokens.push_back("");
    return tokens;
}

// returns 1 if all parsed bodies were already copied, otherwise 0.
i32 IO::Parser::next_body_info(IO::Body& io_body_ref){
    std::string line;
    if (failed)
        throw std::string("Parser failed.");

    if (getline(input_file, line)){
        const auto tokens = tokenize_string(line, '\t');
        if (tokens.size() != 6)
            throw std::string("Body " + std::to_string(parsed_bodies) + " has "
                + std::to_string(tokens.size()) + " tokens.");

        try {
            io_body_ref.id = std::stoi(tokens[0]);
            io_body_ref.mass = std::stod(tokens[1]);
            io_body_ref.x = std::stod(tokens[2]);
            io_body_ref.y = std::stod(tokens[3]);
            io_body_ref.vel_x = std::stod(tokens[4]);
            io_body_ref.vel_y = std::stod(tokens[5]);
        }
        catch (const std::invalid_argument& e){
            failed = true;
            throw std::string("Body " + std::to_string(parsed_bodies) + ": invalid arg.");
        }
        catch (const std::out_of_range& e){
            failed = true;
            throw std::string("Body " + std::to_string(parsed_bodies) + ": out of range.");
        }
    }
    else
        return 1;
    parsed_bodies++;
    return 0;
}
