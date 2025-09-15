#pragma once

#include <fstream>
#include <iostream>
#include <regex>
#include <string>
#include <filesystem>

void updateConfig(const std::string& filename, const std::string& varName, const std::string& newValue);