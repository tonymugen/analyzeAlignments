/*
 * Copyright (c) 2023 Anthony J. Greenberg
 *
 * Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS
 * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
 * IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
 * THE POSSIBILITY OF SUCH DAMAGE.
 */

/// Identify homozygosity runs
/** \file
 * \author Anthony J. Greenberg
 * \copyright Copyright (c) 2023
 * \version 0.1
 *
 * Read a FASTA alignment file and identify low-diversity regions.
 *
 */

#include <fstream>
#include <iostream>

#include "extraFunctions.hpp"
#include "fastaParser.hpp"

int main(int argc, char *argv[]) {
	const std::string cliHelp = "Available command line flags (in any order):\n"
		"  --input-file      file_name (input file name; required).\n"
		"  --window-size     window_size (window size for similarity estimates; required).\n"
		"  --step-size       step_size (step size for similarity estimates; required).\n"
		"  --impute-missing  if set (with no value) replaces missing values with the consensus nucleotide.\n"
		"  --out-file        file_name (output file name; required).\n";
	try {
		std::unordered_map <std::string, std::string> clInfo;
		std::unordered_map <std::string, std::string> stringVariables;
		std::unordered_map <std::string, int> intVariables;
		BayesicSpace::parseCL(argc, argv, clInfo);
		BayesicSpace::extractCLinfo(clInfo, intVariables, stringVariables);
		BayesicSpace::ParseFASTA fastaAlign( stringVariables.at("input-file") );
		if (stringVariables.at("impute-missing") == "set") {
			fastaAlign.imputeMissing();
		}
		size_t windowSize{0};
		if (intVariables.at("window-size") > 0) {
			windowSize = static_cast<size_t>( intVariables.at("window-size") );
		} else {
			throw std::string("ERROR: window size must be > 0");
		}
		size_t stepSize{0};
		if (intVariables.at("window-size") > 0) {
			stepSize = static_cast<size_t>( intVariables.at("step-size") );
		} else {
			throw std::string("ERROR: step size must be > 0");
		}
		auto result{fastaAlign.diversityInWindows(windowSize, stepSize)};
		std::fstream outStream;
		outStream.open(stringVariables.at("out-file"), std::ios::out);
		BayesicSpace::saveDiversityTable(result, outStream);
		outStream.close();
	} catch(std::string &problem) {
		std::cerr << problem << "\n";
		std::cerr << cliHelp;
		return 1;
	}
}
