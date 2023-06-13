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

#include <iterator>
#include <vector>
#include <unordered_map>
#include <utility> // for std::pair
#include <string>
#include <fstream>
#include <algorithm>

#include "fastaParser.hpp"
#include "ssw_cpp.h"

#include <iostream>

using namespace BayesicSpace;

ParseFASTA::ParseFASTA(const std::string &fastaFileName) {
	std::fstream fastaFile;
	std::string fastaLine;
	fastaFile.open(fastaFileName, std::ios::in);
	// get the first line and examine it (skip any empty lines)
	while ( std::getline(fastaFile, fastaLine) && fastaLine.empty() ) {
	}
	if ( fastaLine.empty() || fastaFile.eof() ) {
		throw std::string("ERROR: all lines in ") + fastaFileName + std::string(" are empty in ") +
			std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	if (fastaLine[0] != '>') {
		throw std::string("ERROR: file ") + fastaFileName + std::string(" does not appear to be a FASTA file (no > on the first line) in ") +
			std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	fastaLine.erase(0, 1);                                                                               // erase the ">" at the beginning
	const auto firstNonSpace = fastaLine.find_first_not_of(' ');
	if (firstNonSpace == std::string::npos) {
		throw std::string("ERROR: some non-space characters required in a FASTA header in ") +
			std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	fastaLine.erase(0, firstNonSpace);
	fastaAlignment_.emplace_back(fastaLine, "");
	while ( std::getline(fastaFile, fastaLine) ) {
		if ( fastaLine.empty() ) {
			continue;
		}
		if (fastaLine[0] == '>') {
			fastaLine.erase(0, 1);                                                                       // erase the ">" at the beginning
			const auto locFirstNonSpace = fastaLine.find_first_not_of(' ');
			if (locFirstNonSpace == std::string::npos) {
				throw std::string("ERROR: some non-space characters required in a FASTA header in ") +
					std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
			}
			fastaLine.erase(0, locFirstNonSpace);
			fastaAlignment_.emplace_back(fastaLine, "");
		} else {
			fastaAlignment_.back().second += fastaLine;
		}
	}
	if (fastaAlignment_.size() < 2) {
		throw std::string("ERROR: alignment file ") + fastaFileName + std::string(" must have at least two sequence records in ") +
			std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	const size_t alignmentSize = fastaAlignment_[0].second.size();
	for (const auto &oneElement : fastaAlignment_) {
		if (oneElement.second.size() != alignmentSize) {
			throw std::string("ERROR: all sequences in file ") + fastaFileName + std::string(" must be the same length in ") +
				std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
		}
	}
	fastaFile.close();
	makeConsensus_();
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
		consensus_      = toCopy.consensus_;
	}
	return *this;
}

ParseFASTA& ParseFASTA::operator=(ParseFASTA &&toMove) noexcept {
	if (this != &toMove) {
		fastaAlignment_ = std::move(toMove.fastaAlignment_);
		consensus_      = std::move(toMove.consensus_);
	}
	return *this;
}

std::string ParseFASTA::extractConsensusWindow(const size_t &startIdx, const size_t &windowLength) const {
	std::string window;
	auto first = consensus_.cbegin() + static_cast<std::string::difference_type>(startIdx);
	std::copy_n( first, windowLength, std::back_inserter(window) );
	return window;
}

std::vector< std::pair< size_t, std::vector<uint32_t> > > ParseFASTA::diversityInWindows(const size_t &windowSize, const size_t &stepSize) const {
	std::vector< std::pair< size_t, std::vector<uint32_t> > > result;
	size_t windowStart{0};
	size_t windowEnd{windowSize};
	while ( windowEnd < this->alignmentLength() ) {
		std::unordered_map<std::string, uint32_t> sequenceTable;
		for (const auto &eachSeq : fastaAlignment_) {
			++sequenceTable[eachSeq.second.substr(windowStart, windowSize)];
		}
		std::vector<uint32_t> counts;
		counts.reserve( sequenceTable.size() );
		for (const auto &eachSequence : sequenceTable) {
			counts.push_back(eachSequence.second);
		}
		result.emplace_back( windowStart, std::move(counts) );
		windowStart += stepSize;
		windowEnd   += stepSize;
	}
	return result;
}

std::unordered_map<std::string, uint32_t> ParseFASTA::extractWindow(const size_t &windowStartPosition, const size_t &windowSize) const {
	std::unordered_map<std::string, uint32_t> result;
	for (const auto &eachSeq : fastaAlignment_) {
		++result[eachSeq.second.substr(windowStartPosition, windowSize)];
	}
	return result;
}

AlignmentStatistics ParseFASTA::extractSequence(const std::string &querySequence) const {
	static const int32_t minMaskLen{15};
	int32_t maskLen{static_cast<int32_t>(querySequence.size() / 2)};
	maskLen = maskLen < minMaskLen ? minMaskLen : maskLen;
	StripedSmithWaterman::Aligner aligner;
	StripedSmithWaterman::Filter filter;
	StripedSmithWaterman::Alignment alignment;
	aligner.Align(querySequence.c_str(), consensus_.c_str(), static_cast<int32_t>( consensus_.size() ), filter, &alignment, maskLen);
	if (alignment.ref_begin < 0) {
		throw std::string("ERROR: matching reference start value cannot be negative in ") +
			std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	if (alignment.ref_end < alignment.ref_begin) {
		throw std::string("ERROR: matching reference end must be greater than start in ") +
			std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	if (alignment.query_begin < 0) {
		throw std::string("ERROR: query start value cannot be negative in ") +
			std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	if (alignment.query_end < alignment.query_begin) {
		throw std::string("ERROR: query end must be greater than start in ") +
			std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	AlignmentStatistics result{
		static_cast<size_t>(alignment.ref_begin),
		static_cast<size_t>(alignment.ref_end - alignment.ref_begin),
		static_cast<size_t>(alignment.query_begin),
		static_cast<size_t>(alignment.query_end - alignment.query_begin),
	};
	return result;
}

void ParseFASTA::imputeMissing() {
	const std::string standardNucleotides("AaCcTtGg-");
	for (auto &eachSeq : fastaAlignment_) {
		std::transform(
			eachSeq.second.cbegin(), eachSeq.second.cend(),
			consensus_.cbegin(),
			eachSeq.second.begin(),
			[&standardNucleotides](char nuc1, char nuc2){
				return standardNucleotides.find_first_of(nuc1)== std::string::npos ? nuc2 : nuc1;
			});
	}
}

void ParseFASTA::makeConsensus_() {
	const size_t alignLength = this->alignmentLength();
	const std::string standardNucleotides("AaCcTtGgNn-");
	for (size_t iNuc = 0; iNuc < alignLength; ++iNuc) {
		std::unordered_map<char, uint32_t> nucleotides;
		for (const auto &eachSeq : fastaAlignment_) {
			char curNucleotide{eachSeq.second.at(iNuc)};
			const auto cnPos = standardNucleotides.find_first_of(curNucleotide);
			if (cnPos != std::string::npos) {
				++nucleotides[curNucleotide];
			}
		}
		auto maxCountIt = std::max_element(nucleotides.begin(), nucleotides.end(), 
			[](std::pair<char, uint32_t> count1, std::pair<char, uint32_t> count2){
				return	count1.second < count2.second;
			});
		if (maxCountIt == nullptr) {
			consensus_.push_back('N');
		} else {
			consensus_.push_back(maxCountIt->first);
		}
	}
}
