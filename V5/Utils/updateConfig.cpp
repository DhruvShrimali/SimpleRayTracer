#include <fstream>
#include <iostream>
#include <regex>
#include <string>
#include <filesystem>
#include "updateConfig.hpp"

namespace fs = std::filesystem;

void updateConfig(const std::string& filename, const std::string& varName, const std::string& newValue) {
    fs::path exePath = fs::current_path(); // path where rayTracer.exe is running
    fs::path filePath = exePath / filename;

    if (!fs::exists(filePath)) {
        std::cerr << "File does not exist: " << filePath << "\n";
        return;
    }

    std::ifstream in(filePath);
    if (!in) {
        std::cerr << "Could not open " << filePath << "\n";
        return;
    }

    // std::string output;
    // std::string line;

    // // Build regex dynamically
    // std::string patternStr = "constexpr\\s+[^;\\n]*\\b" + varName + "\\b\\s*=\\s*[^;]*;";
    // std::regex pattern(patternStr);

    // while (std::getline(in, line)) {
    //     if (std::regex_search(line, pattern)) {
    //         // Replace everything after '=' with newValue
    //         std::size_t eqPos = line.find('=');
    //         if (eqPos != std::string::npos) {
    //             line = line.substr(0, eqPos + 1) + " " + newValue + ";";
    //         }
    //     }
    //     output += line + "\n";
    // }
    // in.close();
    
    std::string output;
    std::string line;
    while (std::getline(in, line)) {
        std::size_t pos = line.find(varName);
        if (pos != std::string::npos) {
            std::size_t eqPos = line.find('=', pos);
            if (eqPos != std::string::npos) {
                line = line.substr(0, eqPos + 1) + " " + newValue + ";";
            }
        }
        output += line + "\n";
    }
    in.close();

    std::ofstream out(filename);
    if (!out) {
        std::cerr << "Could not write to " << filename << "\n";
        return;
    }
    out << output;
}


