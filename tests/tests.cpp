/*
 * Copyright (c) 2013 Anthony J. Greenberg
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

/// Tests
/** \file
 * \author Anthony J. Greenberg
 * \copyright Copyright (c) 2023
 * \version 0.1
 *
 * Tests using Catch2.
 *
 */
#include <string>
#include <vector>
#include <utility>
#include <numeric>
#include <algorithm>
#include <unordered_map>
#include <fstream>

#include <iostream>

#include "catch2/catch_test_macros.hpp"
#include "fastaParser.hpp"

TEST_CASE("A FASTA file is properly parsed", "[parser]") { // NOLINT
	const std::string testFASTAfile("../tests/testK.fasta");
	const std::string emptyFASTA("../tests/empty.fasta");
	constexpr size_t trueSeqNum{19};
	constexpr size_t trueAlgnLen{10040};
	BayesicSpace::ParseFASTA testParser(testFASTAfile);
	const auto nSequences  = testParser.sequenceNumber();
	const auto alignLength = testParser.alignmentLength();
	SECTION("Constructor tests") {
		REQUIRE(nSequences == trueSeqNum);
		REQUIRE(alignLength == trueAlgnLen);
		REQUIRE_THROWS( BayesicSpace::ParseFASTA(emptyFASTA) );
	}
	SECTION("FASTA summary tests") {
		constexpr size_t windowStart{600};
		constexpr size_t windowSize{100};
		constexpr size_t stepSize{50};
		constexpr size_t tooBigStart{trueAlgnLen * 2};
		const auto consensus = testParser.extractConsensusWindow(windowStart, windowSize);
		REQUIRE(consensus == "TGTGACTCTGGTAACTAGAGATCCCTCAGACCCTTTTAGTCAGTGTGGAA-AATCTCTAGCAGTGGCGCCCGAACAGGGA-CTTGAAAGCGAAAGTGAAA");
		REQUIRE_THROWS( testParser.extractConsensusWindow(tooBigStart, windowSize) );
		const auto diversity = testParser.diversityInWindows(windowSize, stepSize);
		REQUIRE(diversity.back().first < alignLength);
		std::vector<uint32_t> sumNseq;
		sumNseq.reserve( diversity.size() );
		for (const auto &window : diversity) {
			sumNseq.emplace_back(std::accumulate(window.second.begin(), window.second.end(), uint32_t{0}) );
		}
		const auto minMaxSeqCount = std::minmax_element( sumNseq.begin(), sumNseq.end() );
		// The sum of sequence counts must be the same for all windows and equal to the total number of sequences
		REQUIRE( (*minMaxSeqCount.first) == (*minMaxSeqCount.second) );
		REQUIRE((*minMaxSeqCount.first) == nSequences);
		const auto windowFromPosition = testParser.extractWindow(windowStart, windowSize);
		// unique sequence counts must sum to the total number
		uint32_t sumUnq{0};
		for (const auto &eachSeq : windowFromPosition) {
			sumUnq += eachSeq.second;
		}
		REQUIRE(sumUnq == nSequences);
		REQUIRE_THROWS( testParser.extractWindow(tooBigStart, windowSize) );
		std::fstream fastaQueryFile;
		std::string fastaQueryLine;
		fastaQueryFile.open(std::string("../tests/querySequence.fasta"), std::ios::in);
		std::getline(fastaQueryFile, fastaQueryLine); 
		std::string querySequence;
		while ( std::getline(fastaQueryFile, fastaQueryLine) ) {
			querySequence += fastaQueryLine;
		}
		fastaQueryFile.close();
		const auto queryWindow = testParser.extractSequence(querySequence);
		REQUIRE(queryWindow.referenceStart < alignLength);
		REQUIRE(queryWindow.referenceLength <= alignLength);
	}
}

