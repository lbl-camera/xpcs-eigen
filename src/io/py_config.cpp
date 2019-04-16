#include <boost/python.hpp>
#include <boost/python/object.hpp>
#include <boost/python/dict.hpp>
#include <boost/python/numpy.hpp>

#include <xpcs/configuration.h>

namespace np = boost::python::numpy;
namespace python = boost::python;

template<typename T>
T * getNPData(python::object np_array){
    np::ndarray a = python::extract<np::ndarray>(np_array);
    return reinterpret_cast<T *>(a.get_data());
}

void xpcs::Configuration::from_pydict(python::dict &xpcs) {
    dqmap = getNPData<int>(xpcs.get("dqmap"));
    sqmap = getNPData<int>(xpcs.get("sqmap"));

    frameStart = python::extract<int>(xpcs.get("data_begin"));
    frameEnd = python::extract<int>(xpcs.get("data_end"));
    frameStartTodo = python::extract<int>(xpcs.get("data_begin_todo"));
    frameEndTodo = python::extract<int>(xpcs.get("data_end_todo"));
    delays_per_level_ = python::extract<int>(xpcs.get("delays_per_level"));
    darkFrameStart = python::extract<int>(xpcs.get("dark_begin_todo"));
    darkFrameEnd = python::extract<int>(xpcs.get("dark_end_todo"));

    two2one_window_size_ =
        python::extract<int>(xpcs.get("twotime2onetime_window_size"));
    if (!two2one_window_size_) two2one_window_size_ = 1;

    frame_stride_ = python::extract<long>(xpcs.get("stride_frames"));
    frame_average_ = python::extract<long>(xpcs.get("avg_frames"));

    normalizedByFramesum = python::extract<bool>(xpcs.get("normalize_by_framesum"));
    if (!normalizedByFramesum) normalizedByFramesum = false;

    if (darkFrameStart == darkFrameEnd || darkFrameEnd == 0) {
        darkFrameStart = 0;
        darkFrameEnd = 0;
        darkFrames = 0;
    } else {
        darkFrameStart = darkFrameStart - 1;
        darkFrameEnd = darkFrameEnd - 1;
        darkFrames = darkFrameEnd - darkFrameStart + 1;

        darkThreshold = python::extract<float>(xpcs.get("lld"));
        darkSigma = python::extract<float>(xpcs.get("sigma"));
    }

    m_totalStaticPartitions = 0;
    m_totalDynamicPartitions = 0;

    // detector shit
    python::dict detector = python::extract<python::dict>(xpcs.get("detector"));
    xdim = python::extract<int>(detector.get("x_dimension"));
    ydim = python::extract<int>(detector.get("y_dimension"));
    m_detDpixX = python::extract<float>(detector.get("x_pixel_size"));
    m_detDpixY = python::extract<float>(detector.get("y_pixel_size"));
    m_detAdhupPhot = python::extract<float>(detector.get("adu_per_photon"));
    m_detPreset = python::extract<float>(detector.get("exposure_time"));
    m_detEfficiency = python::extract<float>(detector.get("efficiency"));

    float det_dist = python::extract<float>(detector.get("distance"));
    float flux = python::extract<float>(detector.get("beam_intensity_transmitted"));
    float thickness = python::extract<float>(xpcs.get("thickness"));

    m_normFactor = 1. / m_detEfficiency / m_detAdhupPhot / m_detPreset;
    m_normFactor /= (m_detDpixX / det_dist * m_detDpixY / det_dist);

    m_staticWindow = python::extract<int>(xpcs.get("static_mean_window_size"));
    flatfieldEnabled = python::extract<bool>(xpcs.get("flatfield_enabled"));
    if (flatfieldEnabled) {
        flatfield = getNPData<double>(detector.get("flatfield"));
    } else {
        flatfield = new double[xdim * ydim];
        for (int i = 0; i < xdim * ydim; i++) flatfield[i] = 1.;
    }

    twotime_ = false; // TODO move default values to constructor
    std::string str = python::extract<std::string>(xpcs.get("analysis_type"));
    if (str.compare("Twotime") == 0) {
        twotime_ = true;
        np::ndarray qphi_bins = python::extract<np::ndarray>(xpcs.get("qphi_bin_to_process"));
        int size = qphi_bins.shape(0);
        int * qphibins = reinterpret_cast<int *>(qphi_bins.get_data());
        for (int i = 0; i < size; i++) 
            qphi_bin_to_process_.insert(qphi_bin_to_process_.begin(), qphibins, qphibins + size);
    } 
}
