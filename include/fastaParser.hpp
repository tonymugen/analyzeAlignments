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

/// Class definitions for FASTA alignment parsing
/** \file
 * \author Anthony J. Greenberg
 * \copyright Copyright (c) 2023
 * \version 0.1
 *
 * Class for reading, parsing, and manipulating DNA sequence alignments in FASTA format.
 *
 */

#pragma once

#include <vector>
#include <unordered_map>
#include <utility> // for std::pair
#include <string>
#include <cstdint>

namespace BayesicSpace {
	class ParseFASTA;

	/** \brief FASTA alignment parser
	 *
	 * Reads a FASTA alignment file, separates the sequences and headers, and provides analysis methods.
	 * The data are stored in memory, so users should pay attention to file sizes.
	 *
	 */
	class ParseFASTA {
	public:
		/** \brief Default constructor */
		ParseFASTA() = default;
		/** \brief Constructor from FASTA file 
		 *
		 * Read data from a FASTA file.
		 *
		 * \param[in] fastaFileName input FASTA file name
		 */
		ParseFASTA(const std::string &fastaFileName);
		/** \brief Copy constructor 
		 *
		 * \param[in] toCopy object to copy
		 */
		ParseFASTA(const ParseFASTA &toCopy);
		/** \brief Move constructor 
		 *
		 * \param[in] toMove object to move
		 */
		ParseFASTA(ParseFASTA &&toMove) noexcept;
		/** \brief Copy assignment operator 
		 *
		 * \param[in] toCopy object to copy
		 */
		ParseFASTA& operator=(const ParseFASTA &toCopy);
		/** \brief Move assignment operator 
		 *
		 * \param[in] toMove object to move
		 */
		ParseFASTA& operator=(ParseFASTA &&toMove) noexcept;
		/** \brief Destructor */
		~ParseFASTA() = default;
		/** \brief Number of sequences in alignment 
		 *
		 * \return number of sequences in the alignment
		 */
		size_t sequenceNumber() const noexcept {return fastaAlignment_.size(); };
		/** \brief Alignment length
		 *
		 * \return alignment length
		 */
		size_t alignmentLength() const {return fastaAlignment_.at(0).second.size(); };
		/** \brief Sequence diversity in windows 
		 *
		 * Calculate the number of different sequences in window sliding along a sequence alignment.
		 * Reports the number of times each unique sequence occurs by window position.
		 *
		 * \param[in] windowSize window size in base pairs
		 * \param[in] stepSize window movement steps in base pairs
		 * \return vector of pairs that contain window start positions and unique sequence counts
		 */
		std::vector< std::pair< size_t, std::vector<uint32_t> > > diversityInWindows(const size_t &windowSize, const size_t &stepSize);
		/** \brief Extract an alignment window
		 *
		 * Calculate the number of different sequences in a window.
		 * Reports the number of times each unique sequence occurs in the provided window.
		 *
		 * \param[in] windowStartPosition window start
		 * \param[in] windowSize window size in base pairs
		 * \return map of sequences to the number of times each occurs in the alignment
		 */
		std::unordered_map<std::string, uint32_t> extractWindow(const size_t &windowStartPosition, const size_t &windowSize);
		/** \brief Impute missing values
		 *
		 * Replaces missing (N or other variants, e.g. Y, S, etc.) nucleotides with the consensus value.
		 */
		void imputeMissing();
	private:
		/** \brief Alignment data 
		 *
		 * Each element in the vector is a sequence in the alignment.
		 * The first string in the pair is the FASTA header, the second is the sequence without line breaks.
		 */
		std::vector< std::pair<std::string, std::string> > fastaAlignment_;
	};
}
