#include "fixed.h"
#include "simulation.h"
#include <cstdio>
#include <fstream>
#include <string>

constexpr size_t N = 36, M = 84;

string load_field(const string &filename) {
  ifstream file(filename);
  if (!file.is_open()) {
    throw runtime_error("Failed to open file: " + filename);
  }

  string field;
  string line;
  while (getline(file, line)) {
    field += line + '\n';
  }
  return field;
}

#define FLOAT float
#define DOUBLE double
#define FAST_FIXED(N, K) FastFixed<N, K>
#define FIXED(N, K) Fixed<N, K>

#ifndef TYPES
#define TYPES FLOAT, DOUBLE, FIXED(32, 16), FAST_FIXED(32, 16)
#endif

template <size_t N, size_t M> struct Size {
  static constexpr size_t n = N;
  static constexpr size_t m = M;
};

#define S(N, M) Size<N, M>

#ifndef SIZES
#define SIZES S(36, 84), S(100, 100)
#endif

template <typename T> struct is_fixed : std::false_type {};

template <size_t N, size_t K> struct is_fixed<Fixed<N, K>> : std::true_type {};

template <typename T> struct is_fast_fixed : std::false_type {};

template <size_t N, size_t K>
struct is_fast_fixed<FastFixed<N, K>> : std::true_type {};

template <typename T> static bool is_type(std::string type) {
  if constexpr (std::is_same_v<T, float>) {
    return type == "FLOAT";
  } else if constexpr (std::is_same_v<T, double>) {
    return type == "DOUBLE";
  } else if constexpr (is_fixed<T>::value) {
    if (!type.starts_with("FIXED")) {
      return false;
    }
    int N = 0, K = 0;
    sscanf(type.c_str(), "FIXED(%d,%d)", &N, &K);
    if (T::bits == N && T::frac_bits == K) {
      return true;
    }
  } else if constexpr (is_fast_fixed<T>::value) {
    if (!type.starts_with("FAST_FIXED")) {
      return false;
    }
    int N = 0, K = 0;
    sscanf(type.c_str(), "FAST_FIXED(%d,%d)", &N, &K);
    if (T::bits == N && T::frac_bits == K) {
      return true;
    }
  }
  return false;
}

template <typename T> std::string get_type_name() {
  if (std::is_same_v<T, float>) {
    return "FLOAT";
  } else if (std::is_same_v<T, double>) {
    return "DOUBLE";
  } else if constexpr (is_fixed<T>::value) {
    return "FIXED(" + std::to_string(T::bits) + "," +
           std::to_string(T::frac_bits) + ")";
  } else if constexpr (is_fast_fixed<T>::value) {
    return "FAST_FIXED(" + std::to_string(T::bits) + "," +
           std::to_string(T::frac_bits) + ")";
  } else {
    return "UNKNOWN";
  }
}

template <typename PT, typename VT, typename VFT, size_t N, size_t M>
void run_simulation(std::string field) {
  Simulation<PT, VT, VFT, N, M> sim;
  using State = typename Simulation<PT, VT, VFT, N, M>::State;
  PT rho[256];
  rho[' '] = 0.01;
  rho['.'] = 1000;

  auto state = State(rho, 0.1, field);
  sim.init(state);
  for (size_t i = 0; i < 1'000'000; ++i) {
    sim.step(i);
  }
}

std::pair<size_t, size_t> parse_size(std::string field) {
  size_t N = 0, M = 0;
  for (size_t i = 0; i < field.size(); ++i) {
    if (field[i] == '\n') {
      M = i;
      break;
    }
  }
  N = field.size() / (M + 1);
  return {N, M};
}

template <typename T> static bool is_size(size_t N, size_t M) {
  return T::n == N && T::m == M;
}

template <typename PT, typename VT, typename VFT, typename T, typename... Ts>
bool match_field_size(std::string field) {
  auto [N_, M_] = parse_size(field);
  if (is_size<T>(N_, M_)) {
    cout << "Matched: " << get_type_name<PT>() << " " << get_type_name<VT>()
         << " " << get_type_name<VFT>() << " " << N_ << " " << M_ << endl;
    run_simulation<PT, VT, VFT, T::n, T::m>(field);
    return true;
  } else if constexpr (sizeof...(Ts) > 0) {
    return match_field_size<PT, VT, VFT, Ts...>(field);
  } else {
    return false;
  }
}

template <typename PT, typename VT, typename T, typename... Ts>
bool match_v_flow_type(std::string v_flow_type, std::string field) {
  if (is_type<T>(v_flow_type)) {
    return match_field_size<PT, VT, T, SIZES>(field);
    return true;
  } else if constexpr (sizeof...(Ts) > 0) {
    return match_v_flow_type<PT, VT, Ts...>(v_flow_type, field);
  } else {
    return false;
  }
}

template <typename PT, typename T, typename... Ts>
bool match_v_type(std::string v_type, std::string v_flow_type,
                  std::string field) {
  if (is_type<T>(v_type)) {
    return match_v_flow_type<PT, T, TYPES>(v_flow_type, field);
  } else if constexpr (sizeof...(Ts) > 0) {
    return match_v_type<PT, Ts...>(v_type, v_flow_type, field);
  } else {
    return false;
  }
}

template <typename T, typename... Ts>
bool match_p_type(std::string p_type, std::string v_type,
                  std::string v_flow_type, std::string field) {
  if (is_type<T>(p_type)) {
    return match_v_type<T, TYPES>(v_type, v_flow_type, field);
  } else if constexpr (sizeof...(Ts) > 0) {
    return match_p_type<Ts...>(p_type, v_type, v_flow_type, field);
  } else {
    return false;
  }
}

void run_simulation(std::string p_type, std::string v_type,
                    std::string v_flow_type, std::string field) {
  if (!match_p_type<TYPES>(p_type, v_type, v_flow_type, field)) {
    cout << "No matching types found or size found" << endl;
    return;
  }
}

int main(int argc, char *argv[]) {
  auto parse_arg = [&argc, &argv](std::string arg) {
    for (int i = 1; i < argc; ++i) {
      if (argv[i] == arg) {
        return std::string(argv[i + 1]);
      }
    }
    throw std::runtime_error("Argument not found: " + arg);
  };

  auto field = load_field(parse_arg("--field"));
  auto p_type = parse_arg("--p-type");
  auto v_type = parse_arg("--v-type");
  auto v_flow_type = parse_arg("--v-flow-type");

  run_simulation(p_type, v_type, v_flow_type, field);

  return 0;
}