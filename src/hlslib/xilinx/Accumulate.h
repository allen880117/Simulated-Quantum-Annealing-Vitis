/// @author    Johannes de Fine Licht (definelicht@inf.ethz.ch)
/// @copyright This software is copyrighted under the BSD 3-Clause License.

#pragma once

#include "hlslib/xilinx/DataPack.h"
#include "hlslib/xilinx/Operators.h"
#include "hlslib/xilinx/Simulation.h"
#include "hlslib/xilinx/Stream.h"

// This header includes a module that allows pipelined accumulation of data
// types with a non-zero latency on the operation, such as addition of floating
// point numbers, using a minimal amount of arithmetic units.
// The functionality is split into three functions:
//
//   1) Iterate, which accumulates `latency` different inputs, outputting
//      partial results.
//   2) Feedback, which loops partial results back to Iterate to continue
//      accumulation until done, then emits the result.
//   3) Reduce, which finally collapses the `latency` inputs into a single
//      output. 
//
// If latency is not set high enough, the design will hang, as no new values
// will be available from Feedback when Iterate requires them, and the pipeline
// will stall forever. The latency should be at least the sum of the latencies
// of the Iterate and Feedback modules, plus a few cycles for the overhead of
// the FIFOs between them.
// For example, for single precision floating point addition at 300 MHz, the
// latency of Iterate is 10 and the latency of bounce is 2. A latency of 14 is 
// enough to successfully run accumulation.
//
// Due to a bug where Vivado HLS (as of 2017.1) considers all streams static,
// the three functions cannot be factored into a generic function here, as
// calling it multiple times would result in dataflow errors on the internal
// streams. Instead call them as follows:
//
//   hlslib::Stream<T> toFeedback("fromFeedback");
//   hlslib::Stream<T> fromFeedback("fromFeedback");
//   #pragma HLS STREAM variable=fromFeedback depth=latency
//   Stream<T> toReduce("toReduce");
//   HLSLIB_DATAFLOW_FUNCTION(hlslib::AccumulateIterate, in, fromFeedback,
//                            toFeedback);
//   HLSLIB_DATAFLOW_FUNCTION(hlslib::AccumulateFeedback, toFeedback,
//                            fromFeedback, toReduce);
//   HLSLIB_DATAFLOW_FUNCTION(hlslib::AccumulateReduce(toReduce, out);
//
// For convenience, a function for reducing with a single cycle latency
// operation is also included as AccumulateSimple.

namespace hlslib {

template <typename T, class Operator, int latency>
void AccumulateIterate(Stream<T> &input, Stream<T, latency> &fromFeedback,
                       Stream<T> &toFeedback, int size, int iterations) {
AccumulateIterate_Iterations:
  for (int i = 0; i < iterations; ++i) {
  AcumulateIterate_Size:
    for (int j = 0; j < size / latency; ++j) {
    AccumulateIterate_Latency:
      for (int k = 0; k < latency; ++k) {
        #pragma HLS PIPELINE
        #pragma HLS LOOP_FLATTEN
        const T a = input.Pop();
        T b;
        if (j > 0) {
          b = fromFeedback.ReadOptimistic();
        } else {
          b = Operator::identity();
        }
        const auto result = Operator::Apply(a, b);
        toFeedback.Push(result);
      }
    }
  }
}

template <typename T, int latency>
void AccumulateFeedback(Stream<T> &toFeedback, Stream<T, latency> &fromFeedback,
                        Stream<T> &toReduce, int size, int iterations) {
AccumulateFeedback_Iterations:
  for (int i = 0; i < iterations; ++i) {
  AccumulateFeedback_Size:
    for (int j = 0; j < size / latency; ++j) {
    AccumulateFeedback_Latency:
      for (int k = 0; k < latency; ++k) {
        #pragma HLS PIPELINE
        #pragma HLS LOOP_FLATTEN
        const auto read = toFeedback.Pop();
        if (j < size / latency - 1) {
          // Feedback back
          fromFeedback.Push(read);
        } else {
          // Write to output
          toReduce.Push(read);
        }
      }
    }
  }
}

template <typename T, class Operator, int latency>
void AccumulateReduce(Stream<T> &toReduce, Stream<T> &output, int size,
                      int iterations) {
AccumulateReduce_Iterations:
  for (int i = 0; i < iterations; ++i) {
    #pragma HLS PIPELINE II=latency
    T result(Operator::identity());
  AccumulateReduce_Latency:
    for (int j = 0; j < latency; ++j) {
      const auto read = toReduce.Pop();
      result = Operator::Apply(result, read);
    }
    output.Push(result);
  }
}

/// Trivial implementation for single cycle latency
template <typename T, class Operator>
void AccumulateSimple(Stream<T> &in, Stream<T> &out, int size, int iterations) {
AccumulateSimple_Iterations:
  for (int i = 0; i < iterations; ++i) {
    T acc;
  AccumulateSimple_Size:
    for (int j = 0; j < size; ++j) {
      #pragma HLS LOOP_FLATTEN
      #pragma HLS PIPELINE
      const auto a = in.Pop();
      const auto b = (j > 0) ? acc : T(Operator::identity());
      acc = Operator::Apply(a, b);
      if (j == size - 1) {
        out.Push(acc);
      }
    }
  }
}

} // End namespace hlslib
