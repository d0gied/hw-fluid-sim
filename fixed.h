#pragma once
#include <bits/stdc++.h>
#include <cfloat>
#include <cstddef>
#include <cstdint>
#include <type_traits>

namespace {
template <typename T, typename MulT, size_t K, typename ResT> struct FixedBase {
  constexpr FixedBase(int v) : v(v << K) {}
  constexpr FixedBase(long long v) : v(v << K) {}
  constexpr FixedBase(float f) : v(f * (1 << K)) {}
  constexpr FixedBase(double f) : v(f * (1 << K)) {}
  constexpr FixedBase() : v(0) {}
  template <template <size_t, size_t> class F, size_t N, size_t K2>
  explicit constexpr FixedBase(const F<N, K2> &x) {
    if constexpr (K2 == K) {
      v = x.v;
    } else if constexpr (K2 < K) {
      v = (T)x.v << (K - K2);
    } else {
      v = (T)x.v >> (K2 - K);
    }
  }

  constexpr operator double() const { return v / (double)(1 << K); }
  constexpr operator float() const { return v / (float)(1 << K); }

  T v;

  auto operator<=>(const FixedBase &) const = default;
  bool operator==(const FixedBase &) const = default;

  static constexpr ResT from_raw(T x) {
    ResT ret;
    ret.v = x;
    return ret;
  }
  friend std::ostream &operator<<(std::ostream &out, ResT x) {
    return out << x.v / (double)(1 << K);
  }

  friend ResT abs(ResT x) {
    if (x.v < 0) {
      x.v = -x.v;
    }
    return x;
  }
};
} // namespace

template <template <size_t, size_t> class Fixed, size_t N, size_t K>
Fixed<N, K> operator+(const Fixed<N, K> &a, const Fixed<N, K> &b) {
  return Fixed<N, K>::from_raw(a.v + b.v);
}

template <template <size_t, size_t> class Fixed, size_t N, size_t K>
Fixed<N, K> operator-(const Fixed<N, K> &a, const Fixed<N, K> &b) {
  return Fixed<N, K>::from_raw(a.v - b.v);
}

template <template <typename, typename, size_t, typename> class FixedBase,
          typename T, typename MulT, size_t K, typename ResT>
ResT operator*(const FixedBase<T, MulT, K, ResT> &a,
               const FixedBase<T, MulT, K, ResT> &b) {
  return ResT::from_raw(((MulT)a.v * b.v) >> K);
}

template <template <typename, typename, size_t, typename> class FixedBase,
          typename T, typename MulT, size_t K, typename ResT>
ResT operator/(const FixedBase<T, MulT, K, ResT> &a,
               const FixedBase<T, MulT, K, ResT> &b) {
  return ResT::from_raw(((MulT)a.v << K) / b.v);
}

template <template <size_t, size_t> class Fixed, size_t N, size_t K>
Fixed<N, K> operator-(const Fixed<N, K> &a) {
  return Fixed<N, K>::from_raw(-a.v);
}

template <template <size_t, size_t> class Fixed, size_t N, size_t K>
Fixed<N, K> &operator+=(Fixed<N, K> &a, const Fixed<N, K> &b) {
  return a = a + b;
}

template <template <size_t, size_t> class Fixed, size_t N, size_t K>
Fixed<N, K> &operator-=(Fixed<N, K> &a, const Fixed<N, K> &b) {
  return a = a - b;
}

template <template <size_t, size_t> class Fixed, size_t N, size_t K>
Fixed<N, K> &operator*=(Fixed<N, K> &a, const Fixed<N, K> &b) {
  return a = a * b;
}

template <template <size_t, size_t> class Fixed, size_t N, size_t K>
Fixed<N, K> &operator/=(Fixed<N, K> &a, const Fixed<N, K> &b) {
  return a = a / b;
}

template <size_t N, size_t K> struct Fixed;

template <size_t K>
struct Fixed<16, K> : FixedBase<int16_t, int32_t, K, Fixed<16, K>> {
  static constexpr size_t bits = 16;
  static constexpr size_t frac_bits = K;
  using Fixed16Base = FixedBase<int16_t, int32_t, K, Fixed<16, K>>;
  using Fixed16Base::Fixed16Base;
  constexpr Fixed() : Fixed16Base() {}
  constexpr Fixed(int v) : Fixed16Base(v) {}
  constexpr Fixed(long long v) : Fixed16Base(v) {}
  constexpr Fixed(float f) : Fixed16Base(f) {}
  constexpr Fixed(double f) : Fixed16Base(f) {}
};

template <size_t K>
struct Fixed<32, K> : FixedBase<int32_t, int64_t, K, Fixed<32, K>> {
  static constexpr size_t bits = 32;
  static constexpr size_t frac_bits = K;
  using Fixed32Base = FixedBase<int32_t, int64_t, K, Fixed<32, K>>;
  using Fixed32Base::Fixed32Base;
  constexpr Fixed() : Fixed32Base() {}
  constexpr Fixed(int v) : Fixed32Base(v) {}
  constexpr Fixed(long long v) : Fixed32Base(v) {}
  constexpr Fixed(float f) : Fixed32Base(f) {}
  constexpr Fixed(double f) : Fixed32Base(f) {}
};

template <size_t K>
struct Fixed<64, K> : FixedBase<int64_t, __int128, K, Fixed<64, K>> {
  static constexpr size_t bits = 64;
  static constexpr size_t frac_bits = K;
  using Fixed64Base = FixedBase<int64_t, __int128, K, Fixed<64, K>>;
  using Fixed64Base::Fixed64Base;
  constexpr Fixed() : Fixed64Base() {}
  constexpr Fixed(int v) : Fixed64Base(v) {}
  constexpr Fixed(long long v) : Fixed64Base(v) {}
  constexpr Fixed(float f) : Fixed64Base(f) {}
  constexpr Fixed(double f) : Fixed64Base(f) {}
};

#define FAST_FIXED_TYPE(N)                                                     \
  typename std::conditional<                                                   \
      N <= 8, int_fast8_t,                                                     \
      typename std::conditional<                                               \
          N <= 16, int_fast16_t,                                               \
          typename std::conditional<N <= 32, int_fast32_t,                     \
                                    int_fast64_t>::type>::type>::type

template <size_t N, size_t K>
struct FastFixed : FixedBase<FAST_FIXED_TYPE(N), FAST_FIXED_TYPE(N * 2), K,
                             FastFixed<N, K>> {
  static constexpr size_t bits = N;
  static constexpr size_t frac_bits = K;
  using FastFixedBase =
      FixedBase<FAST_FIXED_TYPE(N), FAST_FIXED_TYPE(N * 2), K, FastFixed<N, K>>;
  using FastFixedBase::FastFixedBase;
  constexpr FastFixed() : FastFixedBase() {}
  constexpr FastFixed(int v) : FastFixedBase(v) {}
  constexpr FastFixed(long long v) : FastFixedBase(v) {}
  constexpr FastFixed(float f) : FastFixedBase(f) {}
  constexpr FastFixed(double f) : FastFixedBase(f) {}
};