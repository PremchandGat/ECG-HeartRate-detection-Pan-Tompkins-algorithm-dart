import 'dart:math';
import 'data.dart';

List<double> convolution(List<double> signal) {
  List<double> sig1 = [...signal];
  List<double> sig2 = [for (int i = 0; i < 20; i++) 0.05];
  List<double> conv = [for (int i = 0; i < (sig1.length - sig2.length); i++) 0];
  for (int l = 0; l < conv.length; l++) {
    for (int i = 0; i < sig2.length; i++) {
      conv[l] += sig1[l - i + sig2.length] * sig2[i];
    }
  }

  List<double> d = [for (int i = 0; i < signal.length; i++) 0];
  for (int l = 0; l < conv.length; l++) {
    d[l + 11] = conv[l];
  }
  return d;
}

double geMax(List signal) {
  double largestGeekValue = signal[0];

  for (var i = 0; i < signal.length; i++) {
    // Checking for largest value in the list
    if (signal[i] > largestGeekValue) {
      largestGeekValue = signal[i];
    }
  }
  return largestGeekValue;
}

List diff(List signal) {
  List out = [];
  if (signal.length < 2) {
    return [];
  }
  for (int i = 0; i < signal.length - 1; i++) {
    out.add(signal[i + 1] - signal[i]);
  }
  return out;
}

List unique(List arr) {
  arr.sort();
  List unique_list = [];

  var last_added;

  for (var element in arr) {
    if (element != last_added) {
      unique_list.add(element);
      last_added = element;
    }
  }
  return unique_list;
}

List divideList(List signal, int divider) {
  List out = [];
  for (int i = 0; i < signal.length; i++) {
    out.add(signal[i] / divider);
  }
  return out;
}

double average(List signal) {
  double summ = 0;
  for (int i = 0; i < signal.length; i++) {
    summ += signal[i];
  }
  return summ / signal.length;
}

class PanTompkinsQRS {
  List<double> bandPassFilter(List<double> signal) {
    List<double> sig = List.from(signal);
    List<double> result = <double>[
      for (double i = 0; i < signal.length; i++) i
    ];

    for (int index = 0; index < signal.length; index++) {
      sig[index] = signal[index];

      if (index >= 1) {
        sig[index] += 2 * sig[index - 1];
      }

      if (index >= 2) {
        sig[index] -= sig[index - 2];
      }

      if (index >= 6) {
        sig[index] -= 2 * signal[index - 6];
      }

      if (index >= 12) {
        sig[index] += signal[index - 12];
      }
    }

    result = List.from(sig);

    for (int index = 0; index < signal.length; index++) {
      result[index] = -1 * sig[index];

      if (index >= 1) {
        result[index] -= result[index - 1];
      }

      if (index >= 16) {
        result[index] += 32 * sig[index - 16];
      }

      if (index >= 32) {
        result[index] += sig[index - 32];
      }
    }

    double maxVal = max(result.reduce(max), -result.reduce(min));
    result = result.map((val) => val / maxVal).toList();

    return result;
  }

  List<double> derivative(List<double> signal, double fs) {
    List<double> result = <double>[
      for (double i = 0; i < signal.length; i++) i
    ];

    for (int index = 0; index < signal.length; index++) {
      result[index] = 0;

      if (index >= 1) {
        result[index] -= 2 * signal[index - 1];
      }

      if (index >= 2) {
        result[index] -= signal[index - 2];
      }

      if (index >= 2 && index <= signal.length - 2) {
        result[index] += 2 * signal[index + 1];
      }

      if (index >= 2 && index <= signal.length - 3) {
        result[index] += signal[index + 2];
      }

      result[index] = (result[index] * fs) / 8;
    }
    return result;
  }

  List<double> squaring(List<double> signal) {
    List<double> result = <double>[
      for (double i = 0; i < signal.length; i++) i
    ];

    for (int index = 0; index < signal.length; index++) {
      result[index] = signal[index] * signal[index];
    }

    return result;
  }

  List<double> movingWindowIntegration(List<double> signal, double fs) {
    List<double> result = <double>[
      for (double i = 0; i < signal.length; i++) i
    ];
    int winSize = (0.150 * fs).round();
    double sum = 0;

    for (int j = 0; j < winSize; j++) {
      sum += signal[j] / winSize;
      result[j] = sum;
    }

    for (int index = winSize; index < signal.length; index++) {
      sum += signal[index] / winSize;
      sum -= signal[index - winSize] / winSize;
      result[index] = sum;
    }

    return result;
  }

  (double, List<int>) solve(List<double> signal, int fs) {
    List<double> inputSignal = List.from(signal);
    List<double> bpass = bandPassFilter(inputSignal);
    List<double> der = derivative(bpass, fs.toDouble());
    List<double> sqr = squaring(der);
    List<double> mwin = movingWindowIntegration(sqr, fs.toDouble());
    List<int> peaks = detectPeaks(
        ecgSingal: signal,
        fs: fs,
        integration_signal: mwin,
        band_pass_signal: bpass);
    double heartRate = (60 * fs) / average(diff(peaks.sublist(1)));
    return (heartRate, peaks);
  }

  List<int> detectPeaks(
      {required List<double> ecgSingal,
      required int fs,
      required List<double> integration_signal,
      required List<double> band_pass_signal}) {
    List<int> possible_peaks = [];
    List<int> signal_peaks = [];
    List<int> r_peaks = [];
    double SPKI = 0;
    double SPKF = 0;
    double NPKI = 0;
    double NPKF = 0;
    List rr_avg_one = [];
    double THRESHOLDI1 = 0;
    double THRESHOLDF1 = 0;
    List<double> rr_avg_two = [];
    double THRESHOLDI2 = 0;
    double THRESHOLDF2 = 0;
    int is_T_found = 0;
    double current_slope = 0;
    double previous_slope = 0;
    int window = (0.15 * fs).round();
    List FM_peaks = [];
    List<double> integration_signal_smooth =
        convolution([...integration_signal]);
    List localDiff = diff([...integration_signal_smooth]);
    double RR_LOW_LIMIT = -999;
    double RR_HIGH_LIMIT = 999;
    double RR_MISSED_LIMIT = 0;

    for (int i = 1; i < localDiff.length; i++) {
      if (i - 1 > 2 * fs && localDiff[i - 1] > 0 && localDiff[i] < 0) {
        FM_peaks.add(i - 1);
      }
    }
    for (int index = 0; index < FM_peaks.length; index++) {
      int current_peak = FM_peaks[index];
      int left_limit = max(current_peak - window, 0);
      int right_limit = min(current_peak + window + 1, band_pass_signal.length);
      int max_index = -1;
      double max_value = -999999;
      for (int i = left_limit; i < right_limit; i++) {
        if (band_pass_signal[i] > max_value) {
          max_value = band_pass_signal[i];
          max_index = i;
        }
      }
      if (max_index != -1) {
        possible_peaks.add(max_index);
      }
      if (index == 0 || index > possible_peaks.length) {
        if (integration_signal[current_peak] >= THRESHOLDI1) {
          SPKI = 0.125 * integration_signal[current_peak] + 0.875 * SPKI;
          if (possible_peaks[index] > THRESHOLDF1) {
            SPKF = 0.125 * band_pass_signal[index] + 0.875 * SPKF;
            signal_peaks.add(possible_peaks[index]);
          } else {
            NPKF = 0.125 * band_pass_signal[index] + 0.875 * NPKF;
          }
        } else if ((integration_signal[current_peak] > THRESHOLDI2 &&
                integration_signal[current_peak] < THRESHOLDI1) ||
            (integration_signal[current_peak] < THRESHOLDI2)) {
          NPKI = 0.125 * integration_signal[current_peak] + 0.875 * NPKI;
          NPKF = 0.125 * band_pass_signal[index] + 0.875 * NPKF;
        }
      } else {
        List RRAVERAGE1 = divideList(
            diff(FM_peaks.sublist(max(0, index - 8), index + 1)), fs);
        double rr_one_mean = average(RRAVERAGE1);
        rr_avg_one.add(rr_one_mean);
        double limit_factor = rr_one_mean;

        if (index >= 8) {
          for (double RR in RRAVERAGE1) {
            if (RR > RR_LOW_LIMIT && RR < RR_HIGH_LIMIT) {
              rr_avg_two.add(RR);
              if (rr_avg_two.length == 9) {
                rr_avg_two.removeAt(0);
                limit_factor = average(rr_avg_two);
              }
            }
          }
          if (rr_avg_two.length == 8 || index < 8) {
            RR_LOW_LIMIT = 0.92 * limit_factor;
            RR_HIGH_LIMIT = 1.16 * limit_factor;
            RR_MISSED_LIMIT = 1.66 * limit_factor;
            RR_MISSED_LIMIT = 1.66 * limit_factor;
          }
          if (rr_avg_one[rr_avg_one.length - 1] < RR_LOW_LIMIT ||
              rr_avg_one[rr_avg_one.length - 1] > RR_MISSED_LIMIT) {
            THRESHOLDI1 = THRESHOLDI1 / 2;
            THRESHOLDF1 = THRESHOLDF1 / 2;
          }

          double curr_rr_interval = RRAVERAGE1[RRAVERAGE1.length - 1];
          int search_back_window = (curr_rr_interval * fs).round();
          if (curr_rr_interval > RR_MISSED_LIMIT) {
            left_limit = current_peak - search_back_window + 1;
            right_limit = current_peak + 1;
            int search_back_max_index = -1;
            max_value = -999999;
            for (int i = left_limit; i < right_limit; i++) {
              if (integration_signal[i] > THRESHOLDI1 &&
                  integration_signal[i] > max_value) {
                max_value = integration_signal[i];
                search_back_max_index = i;
              }
            }
            if (search_back_max_index != -1) {
              SPKI = 0.25 *
                      integration_signal[search_back_max_index < 0
                          ? integration_signal.length + search_back_max_index
                          : search_back_max_index] +
                  0.75 * SPKI;
              THRESHOLDI1 = NPKI + 0.25 * (SPKI - NPKI);
              THRESHOLDI2 = 0.5 * THRESHOLDI1;

              left_limit = search_back_max_index - (0.15 * fs).round();
              right_limit = min(band_pass_signal.length, search_back_max_index);

              int search_back_max_index2 = -1;
              max_value = -999999;

              for (int i = left_limit; i < right_limit; i++) {
                if (band_pass_signal[i] > THRESHOLDF1 &&
                    band_pass_signal[i] > max_value) {
                  max_value = band_pass_signal[i];
                  search_back_max_index2 = i;
                }
              }
              // TODO added if condition
              // if (search_back_max_index2 == -1) {
              //   search_back_max_index2 = band_pass_signal.length - 1;
              // }

              if (band_pass_signal[search_back_max_index2 < 0
                      ? integration_signal.length + search_back_max_index2
                      : search_back_max_index2] >
                  THRESHOLDF2) {
                SPKF = 0.25 *
                        band_pass_signal[search_back_max_index2 < 0
                            ? integration_signal.length + search_back_max_index2
                            : search_back_max_index2] +
                    0.75 * SPKF;
                THRESHOLDF1 = NPKF + 0.25 * (SPKF - NPKF);
                THRESHOLDF2 = 0.5 * THRESHOLDF1;
                signal_peaks.add(search_back_max_index2);
              }
            }
          }
          if (integration_signal[current_peak] >= THRESHOLDI1) {
            if (curr_rr_interval > 0.20 &&
                curr_rr_interval < 0.36 &&
                index > 0) {
              current_slope = geMax(diff(integration_signal.sublist(
                  current_peak - (fs * 0.075).round(), current_peak + 1)));

              previous_slope = geMax(diff(integration_signal.sublist(
                  FM_peaks[index - 1] - (fs * 0.075).round(),
                  FM_peaks[index - 1] + 1)));
              if (current_slope < 0.5 * previous_slope) {
                NPKI = 0.125 * integration_signal[current_peak] + 0.875 * NPKI;
                is_T_found = 1;
              }
            }

            if (is_T_found == 0) {
              SPKI = 0.125 * integration_signal[current_peak] + 0.875 * SPKI;

              if (possible_peaks[index] > THRESHOLDF1) {
                SPKF = 0.125 * band_pass_signal[index] + 0.875 * SPKF;
                signal_peaks.add(possible_peaks[index]);
              } else {
                NPKF = 0.125 * band_pass_signal[index] + 0.875 * NPKF;
              }
            }
          } else if ((integration_signal[current_peak] > THRESHOLDI1 &&
                  integration_signal[current_peak] < THRESHOLDI2) ||
              (integration_signal[current_peak] < THRESHOLDI1)) {
            NPKI = 0.125 * integration_signal[current_peak] + 0.875 * NPKI;
            NPKF = 0.125 * band_pass_signal[index] + 0.875 * NPKF;
          }

          THRESHOLDI1 = NPKI + 0.25 * (SPKI - NPKI);
          THRESHOLDF1 = NPKF + 0.25 * (SPKF - NPKF);
          THRESHOLDI2 = 0.5 * THRESHOLDI1;
          THRESHOLDF2 = 0.5 * THRESHOLDF1;
          is_T_found = 0;
        }
      }
    }
    for (int i in unique(signal_peaks)) {
      int window = (0.2 * fs).round();
      int left_limit = i - window;
      double right_limit = min(i + window + 1, ecgSingal.length.toDouble());
      double max_value = -double.infinity;
      int max_index = -1;
      for (int i = left_limit; i < right_limit; i++) {
        if (ecgSingal[i] > max_value) {
          max_value = ecgSingal[i];
          max_index = i;
        }
      }
      r_peaks.add(max_index);
    }

    return r_peaks;
  }
}

/// vars for High Pass Filter
double hprevFilterd = 0.0;
double hprevUnFiltered = 0.0;
double hprevprevUnfiltered = 0.0;
double hprevprevFilterd = 0.0;

/// Vars for Low Pass Filter
double lprevFilterd = 0.0;
double lprevUnFiltered = 0.0;
double lprevprevUnfiltered = 0.0;
double lprevprevFilterd = 0.0;

/// Low Pass filter
double applyLowPassFilter(double val) {
  double y = 0.2564056711091054 * val +
      0.14992656822522105 * lprevFilterd +
      0.5128113422182107 * lprevUnFiltered -
      0.1755492526616423 * lprevprevFilterd +
      0.2564056711091051 * lprevprevUnfiltered;
  lprevprevFilterd = lprevFilterd;
  lprevFilterd = y;
  lprevprevUnfiltered = lprevUnFiltered;
  lprevUnFiltered = val;
  return y;
}

/// High Pass filter
applyHighPassFilter(double val) {
  double y = 0.9736978852077434 * val +
      1.9467038494842983 * hprevFilterd +
      -1.9473957704154865 * hprevUnFiltered +
      -0.9480876913466759 * hprevprevFilterd +
      0.9736978852077433 * hprevprevUnfiltered;
  hprevprevFilterd = hprevFilterd;
  hprevFilterd = y;
  hprevprevUnfiltered = hprevUnFiltered;
  hprevUnFiltered = val;
  return y;
}

void main(List<String> args) {
  PanTompkinsQRS qrsDetector = PanTompkinsQRS();
  late List<int> peaks;
  late double heartRate;
  List<double> filterdData = [for (double i in data) applyHighPassFilter(i)];
  filterdData = [for (double i in filterdData) applyLowPassFilter(i)];
  (heartRate, peaks) = qrsDetector.solve(filterdData, samplingFreq);
  print("Total Peaks Detected ${peaks.length}");
  print("Heart Rate: " + heartRate.toString() + "BPM");
}
