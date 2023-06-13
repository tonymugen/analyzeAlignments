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
/// Extract unique sequences from an alignment segment
/** \file
 * \author Anthony J. Greenberg
 * \copyright Copyright (c) 2023
 * \version 0.1
 *
 * Read a FASTA alignment file and, extract a segment, and save to a separate file.
 *
 */

#include <fstream>
#include <iostream>
#include <algorithm>
#include <cctype>
#include <string>

#include "extraFunctions.hpp"
#include "fastaParser.hpp"

int main(int argc, char *argv[]) {
	const std::string cliHelp = "Available command line flags (in any order):\n"
		"  --input-file      file_name (input file name; required).\n"
		"  --start-position  start_position (window start position; defaults to 1, first nucleotide).\n"
		"  --window-size     window_size (window size for similarity estimates; required).\n"
		"  --impute-missing  if set (with no value) replaces missing values with the consensus nucleotide.\n"
		"  --query-sequence  a FASTA file with a query sequence to extract a window containing its best match;\n"
		"                    if provided, the --start-position and --window-size flags are ingnored.\n"
		"  --out-format      output file format (FASTA or TAB case-insensitive; defaults to TAB).\n"
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
		size_t startPosition{0};
		if (stringVariables.at("query-sequence") == "unset") {
			if (intVariables.at("window-size") > 0) {
				windowSize = static_cast<size_t>( intVariables.at("window-size") );
			} else {
				throw std::string("ERROR: window size must be > 0");
			}
			if (intVariables.at("start-position") > 0) {
				startPosition = static_cast<size_t>( intVariables.at("start-position") ) - 1;  // make position base-0
			} else {
				throw std::string("ERROR: start position must be greater than 1");
			}
			std::string consensusWindow{fastaAlign.extractConsensusWindow(startPosition, windowSize)};
			// convert to lower case in-place
			std::transform(stringVariables.at("out-format").begin(), stringVariables.at("out-format").end(),
					stringVariables.at("out-format").begin(), [](unsigned char letter){return std::tolower(letter);});
			auto result{fastaAlign.extractWindow(startPosition, windowSize)};
			std::fstream outStream;
			outStream.open(stringVariables.at("out-file"), std::ios::out);
			BayesicSpace::saveUniqueSequences(result, consensusWindow, stringVariables.at("out-format"), outStream);
			outStream.close();
		} else {
			std::fstream fastaQueryFile;
			std::string fastaQueryLine;
			fastaQueryFile.open(stringVariables.at("query-sequence"), std::ios::in);
			std::getline(fastaQueryFile, fastaQueryLine); 
			if (fastaQueryLine[0] != '>') {
				throw std::string("ERROR: file ") + stringVariables.at("query-sequence") + std::string(" does not appear to be a FASTA file (no > on the first line)");
			}
			std::string querySequence;
			while ( std::getline(fastaQueryFile, fastaQueryLine) ) {
				querySequence += fastaQueryLine;
			}
			fastaQueryFile.close();
			BayesicSpace::AlignmentStatistics windowParams{fastaAlign.extractSequence(querySequence)};
			startPosition = windowParams.referenceStart;
			windowSize    = windowParams.referenceLength;
			querySequence = querySequence.substr(windowParams.queryStart, windowParams.queryLength);
			std::string consensusWindow{fastaAlign.extractConsensusWindow(startPosition, windowSize)};
			// convert to lower case in-place
			std::transform(stringVariables.at("out-format").begin(), stringVariables.at("out-format").end(),
					stringVariables.at("out-format").begin(), [](unsigned char letter){return std::tolower(letter);});
			auto result{fastaAlign.extractWindow(startPosition, windowSize)};
			std::fstream outStream;
			outStream.open(stringVariables.at("out-file"), std::ios::out);
			BayesicSpace::saveUniqueSequences(result, consensusWindow, querySequence, stringVariables.at("out-format"), outStream);
			outStream.close();
		}
	} catch(std::string &problem) {
		std::cerr << problem << "\n";
		std::cerr << cliHelp;
		return 1;
	}
}
