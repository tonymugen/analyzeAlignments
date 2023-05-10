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

/// Implementation of FASTA alignment parsing
/** \file
 * \author Anthony J. Greenberg
 * \copyright Copyright (c) 2023
 * \version 0.1
 *
 * Implements the class for reading, parsing, and manipulating DNA sequence alignments in FASTA format.
 *
 */

#include <vector>
#include <utility> // for std::pair
#include <string>
#include <fstream>

#include "fastaParser.hpp"

using namespace BayesicSpace;

ParseFASTA::ParseFASTA(const std::string &fastaFileName) {
	std::fstream fastaFile;
	std::string fastaLine;
	fastaFile.open(fastaFileName, std::ios::in);
	// get the first line and examine it (skip any empty lines)
	while ( std::getline(fastaFile, fastaLine) && !fastaLine.empty() ) {
	}
	if ( fastaLine.empty() ) {
		throw std::string("ERROR: all lines in ") + fastaFileName + std::string(" are empty in ") +
			std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	if (fastaLine[0] != '>') {
		throw std::string("ERROR: file ") + fastaFileName + std::string(" does not appear to be a FASTA file (no > on the first line) in ") +
			std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	fastaLine.erase(0, 1);                                                                       // erase the ">" at the beginning
	const auto firstNonSpace = fastaLine.find_first_not_of(' ');
	if (firstNonSpace == std::string::npos) {
		throw std::string("ERROR: some non-space characters required in a FASTA header in ") +
			std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	fastaLine.erase(0, firstNonSpace);
	fastaAlignment_.emplace_back(fastaLine, std::string());

	fastaFile.close();
}

ParseFASTA::ParseFASTA(const ParseFASTA &toCopy) {
	*this = toCopy;
}

ParseFASTA::ParseFASTA(ParseFASTA &&toMove) noexcept {
	*this = std::move(toMove);
}

ParseFASTA& ParseFASTA::operator=(const ParseFASTA &toCopy) {
	if (this != &toCopy) {
		fastaAlignment_ = toCopy.fastaAlignment_;
	}
	return *this;
}

ParseFASTA& ParseFASTA::operator=(ParseFASTA &&toMove) noexcept {
	if (this != &toMove) {
		fastaAlignment_ = std::move(toMove.fastaAlignment_);
	}
	return *this;
}
