#include <pybind11/pybind11.h>
#include <Las2PowerLine.h>

namespace py = pybind11;

PYBIND11_MODULE(_C, m) {
     m.def("process_points", &process_points);
}