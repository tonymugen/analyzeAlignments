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
 * Definitions of extra utility functions for the FASTA alignment analysis project.
 *
 */

#pragma once

#include <unordered_map>
#include <vector>
#include <utility> // for std::pair
#include <string>
#include <cstdint>
#include <fstream>

namespace BayesicSpace {
	/** \brief Command line parser
	 *
	 * Maps flags to values. Flags assumed to be of the form `--flag-name value`.
	 *
	 * \param[in] argc size of the `argv` array
	 * \param[in] argv command line input array
	 * \param[out] cli map of tags to values
	 */
	void parseCL(int &argc, char **argv, std::unordered_map<std::string, std::string> &cli);
	/** \brief Extract parameters from parsed command line interface flags
	 *
	 * Extracts needed variable values, indexed by `std::string` encoded variable names.
	 *
	 * \param[in] parsedCLI flag values parsed from the command line
	 * \param[out] intVariables indexed `int` variables for use by `main()`
	 * \param[out] stringVariables indexed `std::string` variables for use by `main()`
	 */
	void extractCLinfo(const std::unordered_map<std::string, std::string> &parsedCLI, std::unordered_map<std::string, int> &intVariables, std::unordered_map<std::string, std::string> &stringVariables);
	/** \brief Save the diversity table 
	 *
	 * Save the diversity table. The output file will have two columns: 
	 *     (1) window start position (repeated for every unique sequence).
	 *     (2) number of unique sequence occurrences.
	 * 
	 * \param[in] diversityTable the diversity table data
	 * \param[in] outFile output file stream
	 */
	void saveDiversityTable(const std::vector< std::pair< size_t, std::vector<uint32_t> > > &diversityTable, std::fstream &outFile);
}
