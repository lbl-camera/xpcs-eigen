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
#ifndef CONFIGURATION_
#define CONFIGURATION_

#include <string.h>
#include <map>
#include <vector>

#if (HAVE_HDF5)
#include <hdf5.h>
#endif

#if (HAVE_BOOST_PYTHON)
#include <boost/python.hpp>
#include <boost/python/dict.hpp>
#endif

#include "benchmark.h"

namespace xpcs  {

class Configuration  {

public:    
  Configuration();

  ~Configuration();

  int* getDQMap() const { return dqmap; }
  int* getSQMap() const { return sqmap; }
  int *PixelsPerStaticBin() const { return pixels_per_bin; }

  int getFrameWidth() const { return xdim; }
  int getFrameHeight() const { return ydim; }

  // Starting index of frames that we need to process. 
  int getFrameStartTodo() const { return frameStartTodo - 1; }
  int getFrameEndTodo() const { return frameEndTodo - 1; }
  int getFrameStart() const { return frameStart - 1; }
  int getFrameEnd() const { return frameEnd - 1; }
  int getFrameTodoCount() const;

  // Total number of frames without any stride or average. 
  int getRealFrameTodoCount() const { return frameEndTodo - frameStartTodo + 1; }
  int getFrameCount() const { return frameEnd - frameStart + 1; }
  int getDarkFrameStart() const { return darkFrameStart; }
  int getDarkFrameEnd() const { return darkFrameEnd; }
  int getDarkFrames() const { return darkFrames; }
  float getDarkThreshold() const { return darkThreshold; }
  float getDarkSigma() const { return darkSigma; }

  std::map<int, std::map<int, std::vector<int>> > getBinMaps() const { return m_mapping; }

  int getTotalStaticPartitions() const { return m_totalStaticPartitions; }
  int getTotalDynamicPartitions() const { return m_totalDynamicPartitions; }
  int getStaticWindowSize() const { return m_staticWindow; }
  int DelaysPerLevel() const { return delays_per_level_; }

  int FrameStride() const { return frame_stride_; }
  int FrameAverage() const { return frame_average_; }

  int Two2OneWindowSize() const { return two2one_window_size_; }
  
  std::string getFilename() const { return m_filename; }
  std::string getIMMFilePath() const { return m_immFile; }
  std::string OutputPath() const { return output_path_; }
  void setIMMFilePath(std::string& path) { m_immFile = path; }

  short* getPixelMask() const { return m_validPixelMask; }
  int* getSbinMask() const { return m_sbin; }
  const std::vector<int>& TwoTimeQMask() const { return qphi_bin_to_process_; }

  float getDetDpixX() const { return m_detDpixX; }
  float getDetDpixY() const { return m_detDpixY; }
  float getDetAdhuPhot() const { return m_detAdhupPhot; }
  float getDetPreset() const { return m_detPreset; }
  float getDetEfficiency() const { return m_detEfficiency; }
  float getNormFactor() const { return m_normFactor; }


  bool getIsFlatFieldEnabled() const { return flatfieldEnabled; }
  double* getFlatField() const { return flatfield; }

  bool getIsCompressionEnabled() const { return compression; }
  bool IsNormalizedByFramesum() const { return normalizedByFramesum; }
  bool IsTwoTime() const { return twotime_; }

#if (HAVE_HDF5)
  void from_hdf5(const std::string &, const std::string &);
#endif

#if (HAVE_BOOST_PYTHON)
  void from_pydict(boost::python::dict &);
#endif

private:


  std::string getString(const std::string &path);

  float getFloat(const std::string &path);

  int getInteger(const std::string &path);

  long getLong(const std::string &path);

  int* Dim2DTable(const std::string &path);

  int* get2DTable(const std::string &path);
  
  long* get2DTableL(const std::string &path);

  double* get2DTableD(const std::string &path);

  void BuildQMap();
  
  double* flatfield;

  // Valid pixel mask - mark an entry as 1 in an array if the pixel is contained in any of the bins.
  short* m_validPixelMask;
  int* m_sbin;

  int m_totalStaticPartitions;
  int m_totalDynamicPartitions;

  // Map of dynamic bins to static bin to pixels.
  std::map<int, std::map<int, std::vector<int> >> m_mapping;

  int xdim;
  int ydim;
  int frameStart;
  int frameEnd;
  int frameStartTodo;
  int frameEndTodo;
  int darkFrameStart;
  int darkFrameEnd;
  int darkFrames;
  int m_staticWindow;
  int delays_per_level_;
  int two2one_window_size_;

  int *pixels_per_bin;
  int *dqmap;
  int *sqmap;
  
  int frame_stride_;
  int frame_average_;

  float m_detDpixX;
  float m_detDpixY;
  float m_detAdhupPhot;
  float m_detPreset;
  float m_detEfficiency;
  float m_normFactor;
  float darkThreshold;
  float darkSigma;

#if HAVE_HDF5
  hid_t file_id;
#endif // HAVE_HDF5

  // Flags for checkig if certain fields are enabled. 
  bool compression;
  bool flatfieldEnabled;
  bool normalizedByFramesum;
  bool twotime_;

  std::string m_filename;
  std::string m_immFile;
  std::string output_path_;

  std::vector<int> qphi_bin_to_process_;

};

} // namespace xpcs

#endif
