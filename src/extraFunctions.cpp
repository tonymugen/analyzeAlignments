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

/// Extra functions
/** \file
 * \author Anthony J. Greenberg
 * \copyright Copyright (c) 2023
 * \version 0.1
 *
 * Implementation of extra utility functions for the FASTA alignment analysis project.
 *
 */

#include <array>
#include <algorithm>
#include <iterator>

#include "extraFunctions.hpp"

using namespace BayesicSpace;

void BayesicSpace::parseCL(int &argc, char **argv, std::unordered_map<std::string, std::string> &cli) {
	// set to true after encountering a flag token (the characters after the dash)
	bool val = false;
	// store the token value here
	std::string curFlag;

	for (int iArg = 1; iArg < argc; iArg++) {
		const char *pchar = argv[iArg];
		if ( (pchar[0] == '-') && (pchar[1] == '-') ) { // encountered the double dash, look for the token after it
			if (val) { // A previous flag had no value
				cli[curFlag] = "set";
			}
			// what follows the dash?
			val     = true;
			curFlag = pchar + 2;
		} else {
			if (val) {
				val          = false;
				cli[curFlag] = pchar;
			}
		}
	}
}

void BayesicSpace::extractCLinfo(const std::unordered_map<std::string, std::string> &parsedCLI, std::unordered_map<std::string, int> &intVariables, std::unordered_map<std::string, std::string> &stringVariables) {
	intVariables.clear();
	stringVariables.clear();
	const std::array<std::string, 2> requiredStringVariables{"input-file", "out-file"};
	const std::array<std::string, 4> optionalStringVariables{"impute-missing", "out-format", "query-sequence", "sorted"};
	const std::array<std::string, 3> optionalIntVariables{"start-position", "window-size", "step-size"};
	const std::unordered_map<std::string, std::string> defaultStringValues{ {"impute-missing", "unset"}, {"out-format", "tab"}, {"query-sequence", "unset"}, {"sorted", "unset"} };
	const std::unordered_map<std::string, int> defaultIntValues{ {"start-position", 1}, {"window-size", 100}, {"step-size", 10} };

	if ( parsedCLI.empty() ) {
		throw std::string("No command line flags specified;");
	}
	for (const auto &eachFlag : optionalIntVariables) {
		try {
			intVariables[eachFlag] = stoi( parsedCLI.at(eachFlag) );
		} catch(const std::exception &problem) {
			intVariables[eachFlag] = defaultIntValues.at(eachFlag);
		}
	}
	for (const auto &eachFlag : requiredStringVariables) {
		try {
			stringVariables[eachFlag] = parsedCLI.at(eachFlag);
		} catch(const std::exception &problem) {
			throw std::string("ERROR: ") + eachFlag + std::string(" specification is required");
		}
	}
	for (const auto &eachFlag : optionalStringVariables) {
		try {
			stringVariables[eachFlag] = parsedCLI.at(eachFlag);
		} catch(const std::exception &problem) {
			stringVariables[eachFlag] = defaultStringValues.at(eachFlag);
		}
	}
}

void BayesicSpace::saveDiversityTable(const std::vector< std::pair< size_t, std::vector<uint32_t> > > &diversityTable, std::fstream &outFile) {
	outFile << "position\tcount\n";
	for (const auto &eachWindow : diversityTable) {
		for (const auto &count : eachWindow.second) {
			outFile << eachWindow.first + 1 << "\t" << count << "\n";
		}
	}
}

void BayesicSpace::saveUniqueSequences(const std::unordered_map<std::string, uint32_t> &uniqueSequences, const std::string &consensus, const std::string &fileType, std::fstream &outFile) {
	if (fileType == "fasta") {
		uint32_t seqIdx{1};
		outFile << "> Consensus\n";
		outFile << consensus << "\n";
		for (const auto &eachSeq : uniqueSequences) {
			std::string diffs;
			std::transform(
				eachSeq.first.cbegin(), eachSeq.first.cend(),
				consensus.cbegin(), 
				std::back_inserter(diffs), 
				[](unsigned char nuc1, unsigned char nuc2){return std::toupper(nuc1) == std::toupper(nuc2) ? '.' : std::toupper(nuc1);});
			outFile << "> Sequence " << seqIdx << ": " << eachSeq.second << "\n";
			outFile << diffs << "\n";
			++seqIdx;
		}
	} else if (fileType == "tab") {
		outFile << consensus << "\t" << "C\n";
		for (const auto &eachSeq : uniqueSequences) {
			std::string diffs;
			std::transform(
				eachSeq.first.cbegin(), eachSeq.first.cend(),
				consensus.cbegin(), 
				std::back_inserter(diffs), 
				[](unsigned char nuc1, unsigned char nuc2){return std::toupper(nuc1) == std::toupper(nuc2) ? '.' : std::toupper(nuc1);});
			outFile << diffs << "\t" << eachSeq.second << "\n";
		}
	} else {
		throw std::string("ERROR: output file format must be fasta or tab in ") +
			std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
}

void BayesicSpace::saveUniqueSequences(const std::vector< std::pair<std::string, uint32_t> > &uniqueSequences, const std::string &consensus, const std::string &fileType, std::fstream &outFile) {
	if (fileType == "fasta") {
		uint32_t seqIdx{1};
		outFile << "> Consensus\n";
		outFile << consensus << "\n";
		for (const auto &eachSeq : uniqueSequences) {
			std::string diffs;
			std::transform(
				eachSeq.first.cbegin(), eachSeq.first.cend(),
				consensus.cbegin(), 
				std::back_inserter(diffs), 
				[](unsigned char nuc1, unsigned char nuc2){return std::toupper(nuc1) == std::toupper(nuc2) ? '.' : std::toupper(nuc1);});
			outFile << "> Sequence " << seqIdx << ": " << eachSeq.second << "\n";
			outFile << diffs << "\n";
			++seqIdx;
		}
	} else if (fileType == "tab") {
		outFile << consensus << "\t" << "C\n";
		for (const auto &eachSeq : uniqueSequences) {
			std::string diffs;
			std::transform(
				eachSeq.first.cbegin(), eachSeq.first.cend(),
				consensus.cbegin(), 
				std::back_inserter(diffs), 
				[](unsigned char nuc1, unsigned char nuc2){return std::toupper(nuc1) == std::toupper(nuc2) ? '.' : std::toupper(nuc1);});
			outFile << diffs << "\t" << eachSeq.second << "\n";
		}
	} else {
		throw std::string("ERROR: output file format must be fasta or tab in ") +
			std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
}

void BayesicSpace::saveUniqueSequences(const std::unordered_map<std::string, uint32_t> &uniqueSequences, const std::string &consensus,
								const AlignmentStatistics &alignStats, const std::string &query,
								const std::string &fileType, std::fstream &outFile) {
	if (fileType == "fasta") {
		uint32_t seqIdx{1};
		outFile << "> Query \n";
		outFile << query << "\n";
		outFile << "> Consensus; start: " << alignStats.referenceStart << "; length: " << alignStats.referenceLength << "\n";
		outFile << consensus << "\n";
		for (const auto &eachSeq : uniqueSequences) {
			std::string diffs;
			std::transform(
				eachSeq.first.cbegin(), eachSeq.first.cend(),
				consensus.cbegin(), 
				std::back_inserter(diffs), 
				[](unsigned char nuc1, unsigned char nuc2){return std::toupper(nuc1) == std::toupper(nuc2) ? '.' : std::toupper(nuc1);});
			outFile << "> Sequence " << seqIdx << ": " << eachSeq.second << "\n";
			outFile << diffs << "\n";
			++seqIdx;
		}
	} else if (fileType == "tab") {
		outFile << query     << "\t" << "Q\n";
		outFile << consensus << "\t" << "C|" << alignStats.referenceStart << "|" << alignStats.referenceLength <<"\n";
		for (const auto &eachSeq : uniqueSequences) {
			std::string diffs;
			std::transform(
				eachSeq.first.cbegin(), eachSeq.first.cend(),
				consensus.cbegin(), 
				std::back_inserter(diffs), 
				[](unsigned char nuc1, unsigned char nuc2){return std::toupper(nuc1) == std::toupper(nuc2) ? '.' : std::toupper(nuc1);});
			outFile << diffs << "\t" << eachSeq.second << "\n";
		}
	} else {
		throw std::string("ERROR: output file format must be fasta or tab in ") +
			std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
}

void BayesicSpace::saveUniqueSequences(const std::vector< std::pair<std::string, uint32_t> > &uniqueSequences, const std::string &consensus,
								const AlignmentStatistics &alignStats, const std::string &query,
								const std::string &fileType, std::fstream &outFile) {
	if (fileType == "fasta") {
		uint32_t seqIdx{1};
		outFile << "> Query \n";
		outFile << query << "\n";
		outFile << "> Consensus; start: " << alignStats.referenceStart << "; length: " << alignStats.referenceLength << "\n";
		outFile << consensus << "\n";
		for (const auto &eachSeq : uniqueSequences) {
			std::string diffs;
			std::transform(
				eachSeq.first.cbegin(), eachSeq.first.cend(),
				consensus.cbegin(), 
				std::back_inserter(diffs), 
				[](unsigned char nuc1, unsigned char nuc2){return std::toupper(nuc1) == std::toupper(nuc2) ? '.' : std::toupper(nuc1);});
			outFile << "> Sequence " << seqIdx << ": " << eachSeq.second << "\n";
			outFile << diffs << "\n";
			++seqIdx;
		}
	} else if (fileType == "tab") {
		outFile << query     << "\t" << "Q\n";
		outFile << consensus << "\t" << "C|" << alignStats.referenceStart << "|" << alignStats.referenceLength <<"\n";
		for (const auto &eachSeq : uniqueSequences) {
			std::string diffs;
			std::transform(
				eachSeq.first.cbegin(), eachSeq.first.cend(),
				consensus.cbegin(), 
				std::back_inserter(diffs), 
				[](unsigned char nuc1, unsigned char nuc2){return std::toupper(nuc1) == std::toupper(nuc2) ? '.' : std::toupper(nuc1);});
			outFile << diffs << "\t" << eachSeq.second << "\n";
		}
	} else {
		throw std::string("ERROR: output file format must be fasta or tab in ") +
			std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
}
