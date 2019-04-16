/**

Copyright (c) 2016, UChicago Argonne, LLC. All rights reserved.

Copyright 2016. UChicago Argonne, LLC. This software was produced
under U.S. Government contract DE-AC02-06CH11357 for Argonne National
Laboratory (ANL), which is operated by UChicago Argonne, LLC for the
U.S. Department of Energy. The U.S. Government has rights to use,
reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR
UChicago Argonne, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR a
ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is
modified to produce derivative works, such modified software should
be clearly marked, so as not to confuse it with the version available
from ANL.

Additionally, redistribution and use in source and binary forms, with
or without modification, are permitted provided that the following
conditions are met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in
      the documentation and/or other materials provided with the
      distribution.

    * Neither the name of UChicago Argonne, LLC, Argonne National
      Laboratory, ANL, the U.S. Government, nor the names of its
      contributors may be used to endorse or promote products derived
      from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY UChicago Argonne, LLC AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL UChicago
Argonne, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.

**/

#include "configuration.h"

#include <iostream>
#include <set>

#include "benchmark.h"
#include "hdf5.h"

namespace xpcs {

Configuration::Configuration() : dqmap(NULL), sqmap(NULL), flatfield(NULL) {}

Configuration::~Configuration() {
    // TODO delete tables.
}

void Configuration::BuildQMap() {
    this->m_validPixelMask = new short[this->xdim * this->ydim];
    m_sbin = new int[this->xdim * this->ydim];

    std::map<int, std::set<int>> sbin_to_qbin;

    // Build mapping from q-bin to s-bins and from s-bins to pixels.
    for (int i = 0; i < (xdim * ydim); i++) {
        if (dqmap[i] < 1 || sqmap[i] < 1) continue;
        // Mark this pixel to be part of a valid mask.
        m_validPixelMask[i] = 1;
        m_sbin[i] = sqmap[i];

        // {qbins} - > {sbin -> [pixels]}
        std::map<int, std::map<int, std::vector<int>>>::iterator it =
            m_mapping.find(dqmap[i]);

        // Highest q-map value is equal to total number of dynamic partitions.
        if (this->m_totalDynamicPartitions < dqmap[i])
            this->m_totalDynamicPartitions = dqmap[i];

        if (it != m_mapping.end()) {
            std::map<int, std::vector<int>>& v = it->second;
            std::map<int, std::vector<int>>::iterator it2 = v.find(sqmap[i]);

            if (this->m_totalStaticPartitions < sqmap[i])
                this->m_totalStaticPartitions = sqmap[i];

            if (it2 != v.end()) {
                std::vector<int>& v2 = it2->second;
                v2.push_back(i);
            } else {
                std::vector<int> data;
                data.push_back(i);
                v[sqmap[i]] = data;
            }
        } else {
            std::map<int, std::vector<int>> mapping;
            std::vector<int> data;
            data.push_back(i);
            mapping[sqmap[i]] = data;
            m_mapping[dqmap[i]] = mapping;
        }

        // Save reverse mapping from sbin to qbin for removing duplicates next.
        std::map<int, std::set<int>>::iterator sbin_it =
            sbin_to_qbin.find(sqmap[i]);
        if (sbin_it != sbin_to_qbin.end()) {
            std::set<int>& avec = sbin_it->second;
            avec.insert(dqmap[i]);
        } else {
            std::set<int> avec;
            avec.insert(dqmap[i]);
            sbin_to_qbin[sqmap[i]] = avec;
        }
    }

    // remove duplicates.
    for (auto it = sbin_to_qbin.begin(); it != sbin_to_qbin.end(); it++) {
        int sbin = it->first;
        std::set<int> qbins = it->second;

        if (qbins.size() > 1) {
            // duplicate sbin -> qbin mapping
            // int dups[qbins.size() - 1] = {0};
            int max_qbin = 0;
            int max_pixels = 0;

            for (auto qid = qbins.begin(); qid != qbins.end(); qid++) {
                int q = *qid;
                auto m1 = m_mapping.find(q)->second;
                std::vector<int> m2 = m1.find(sbin)->second;
                int pixels = m2.size();

                if (pixels > max_pixels) {
                    max_pixels = pixels;
                    max_qbin = q;
                }
            }

            auto m1 = m_mapping.find(max_qbin)->second;
            std::vector<int>& dest_sbin = m1.find(sbin)->second;

            for (auto qid = qbins.begin(); qid != qbins.end(); qid++) {
                int q = *qid;
                if (q == max_qbin) continue;

                std::map<int, std::vector<int>>& m1 = m_mapping.find(q)->second;
                std::vector<int>& src_sbin = m1.find(sbin)->second;

                for (std::vector<int>::iterator qbin_it = src_sbin.begin();
                     qbin_it != src_sbin.end(); qbin_it++) {
                    int p = *qbin_it;
                    dest_sbin.push_back(p);
                }

                m1.erase(m1.find(sbin));
            }
        }
    }

    pixels_per_bin = new int[m_totalStaticPartitions];
    for (int i = 0; i < m_totalStaticPartitions; i++) { pixels_per_bin[i] = 0; }

    for (int i = 0; i < (xdim * ydim); i++) {
        if (sqmap[i] < 1 || dqmap[i] < 1) continue;

        pixels_per_bin[sqmap[i] - 1]++;
    }
}

int Configuration::getFrameTodoCount() const {
    int frames = frameEndTodo - frameStartTodo + 1;

    int steps = frame_stride_ * frame_average_;

    return frames / steps;
}

}  // namespace xpcs
