#pragma once
#include <bits/stdc++.h>
#include <cstddef>

using namespace std;

constexpr std::array<pair<int, int>, 4> deltas{
    {{-1, 0}, {1, 0}, {0, -1}, {0, 1}}};

mt19937 rnd(1337);

template <typename PT, typename VT, typename VFT, size_t N, size_t M>
struct Simulation {
  template <typename T> struct VectorField {
    array<T, deltas.size()> v[N][M];
    T &add(int x, int y, int dx, int dy, T dv) {
      return get(x, y, dx, dy) += dv;
    }

    T &get(int x, int y, int dx, int dy) {
      size_t i = ranges::find(deltas, pair(dx, dy)) - deltas.begin();
      assert(i < deltas.size());
      return v[x][y][i];
    }
  };

  struct State {
    State(PT rho[256], VT g, std::string field) {
      memcpy(this->rho, rho, sizeof(this->rho));
      for (size_t i = 0; i < N; ++i) {
        memcpy(this->field[i], field.data() + i * (M + 1), M);
      }
      this->g = g;

      for (size_t x = 0; x < N; ++x) {
        for (size_t y = 0; y < M; ++y) {
          if (this->field[x][y] == '#')
            continue;
          for (auto [dx, dy] : deltas) {
            dirs[x][y] += (this->field[x + dx][y + dy] != '#');
          }
        }
      }
    }
    PT rho[256];
    PT p[N][M]{}, old_p[N][M];
    VT g;
    char field[N][M + 1];
    VectorField<VT> velocity{};
    VectorField<VFT> velocity_flow{};
    int last_use[N][M]{};
    int UT = 0;
    int dirs[N][M]{};
  };

  PT rho[256];
  PT p[N][M]{}, old_p[N][M];
  VT g;
  char field[N][M + 1];
  VectorField<VT> velocity{};
  VectorField<VFT> velocity_flow{};
  int last_use[N][M]{};
  int UT = 0;

  template <typename T> VT to_vt(T x) {
    if constexpr (is_same_v<T, VT>) {
      return x;
    } else {
      return VT(x);
    }
  }

  template <typename T> VFT to_vft(T x) {
    if constexpr (is_same_v<T, VFT>) {
      return x;
    } else {
      return VFT(x);
    }
  }

  template <typename T> PT to_pt(T x) {
    if constexpr (is_same_v<T, PT>) {
      return x;
    } else {
      return PT(x);
    }
  }

  tuple<VT, bool, pair<int, int>> propagate_flow(int x, int y, VT lim) {
    last_use[x][y] = UT - 1;
    VT ret = 0;
    for (auto [dx, dy] : deltas) {
      int nx = x + dx, ny = y + dy;
      if (field[nx][ny] != '#' && last_use[nx][ny] < UT) {
        auto cap = velocity.get(x, y, dx, dy);
        auto flow = velocity_flow.get(x, y, dx, dy);
        if (flow == to_vft(cap)) {
          continue;
        }
        // assert(v >= velocity_flow.get(x, y, dx, dy));
        auto vp = min(lim, cap - to_vt(flow));
        if (last_use[nx][ny] == UT - 1) {
          velocity_flow.add(x, y, dx, dy, to_vft(vp));
          last_use[x][y] = UT;
          // cerr << x << " " << y << " -> " << nx << " " << ny << " " << vp
          // << " / " << lim << "\n";
          return {vp, 1, {nx, ny}};
        }
        auto [t, prop, end] = propagate_flow(nx, ny, vp);
        ret += t;
        if (prop) {
          velocity_flow.add(x, y, dx, dy, to_vft(t));
          last_use[x][y] = UT;
          // cerr << x << " " << y << " -> " << nx << " " << ny << " " << t <<
          // " / " << lim << "\n";
          return {t, prop && end != pair(x, y), end};
        }
      }
    }
    last_use[x][y] = UT;
    return {ret, 0, {0, 0}};
  }

  VT random01() {
    return VT(static_cast<double>((rnd() & ((1 << 16) - 1))) / (1 << 16));
  }

  void propagate_stop(int x, int y, bool force = false) {
    if (!force) {
      bool stop = true;
      for (auto [dx, dy] : deltas) {
        int nx = x + dx, ny = y + dy;
        if (field[nx][ny] != '#' && last_use[nx][ny] < UT - 1 &&
            velocity.get(x, y, dx, dy) > VT(0)) {
          stop = false;
          break;
        }
      }
      if (!stop) {
        return;
      }
    }
    last_use[x][y] = UT;
    for (auto [dx, dy] : deltas) {
      int nx = x + dx, ny = y + dy;
      if (field[nx][ny] == '#' || last_use[nx][ny] == UT ||
          velocity.get(x, y, dx, dy) > VT(0)) {
        continue;
      }
      propagate_stop(nx, ny);
    }
  }

  VT move_prob(int x, int y) {
    VT sum = 0;
    for (size_t i = 0; i < deltas.size(); ++i) {
      auto [dx, dy] = deltas[i];
      int nx = x + dx, ny = y + dy;
      if (field[nx][ny] == '#' || last_use[nx][ny] == UT) {
        continue;
      }
      auto v = velocity.get(x, y, dx, dy);
      if (v < VT(0)) {
        continue;
      }
      sum += v;
    }
    return sum;
  }

  struct ParticleParams {
    ParticleParams(Simulation &sim) : sim(sim) {}
    char type;
    PT cur_p;
    array<VT, deltas.size()> v;
    Simulation &sim;

    void swap_with(int x, int y) {
      swap(sim.field[x][y], type);
      swap(sim.p[x][y], cur_p);
      swap(sim.velocity.v[x][y], v);
    }
  };

  bool propagate_move(int x, int y, bool is_first) {
    last_use[x][y] = UT - is_first;
    bool ret = false;
    int nx = -1, ny = -1;
    do {
      std::array<VT, deltas.size()> tres;
      VT sum = 0;
      for (size_t i = 0; i < deltas.size(); ++i) {
        auto [dx, dy] = deltas[i];
        int nx = x + dx, ny = y + dy;
        if (field[nx][ny] == '#' || last_use[nx][ny] == UT) {
          tres[i] = sum;
          continue;
        }
        auto v = velocity.get(x, y, dx, dy);
        if (v < VT(0)) {
          tres[i] = sum;
          continue;
        }
        sum += v;
        tres[i] = sum;
      }

      if (sum == VT(0)) {
        break;
      }

      VT p = random01() * sum;
      size_t d = std::ranges::upper_bound(tres, p) - tres.begin();

      auto [dx, dy] = deltas[d];
      nx = x + dx;
      ny = y + dy;
      assert(velocity.get(x, y, dx, dy) > VT(0) && field[nx][ny] != '#' &&
             last_use[nx][ny] < UT);

      ret = (last_use[nx][ny] == UT - 1 || propagate_move(nx, ny, false));
    } while (!ret);
    last_use[x][y] = UT;
    for (size_t i = 0; i < deltas.size(); ++i) {
      auto [dx, dy] = deltas[i];
      int nx = x + dx, ny = y + dy;
      if (field[nx][ny] != '#' && last_use[nx][ny] < UT - 1 &&
          velocity.get(x, y, dx, dy) < VT(0)) {
        propagate_stop(nx, ny);
      }
    }
    if (ret) {
      if (!is_first) {
        ParticleParams pp(*this);
        pp.swap_with(x, y);
        pp.swap_with(nx, ny);
        pp.swap_with(x, y);
      }
    }
    return ret;
  }

  int dirs[N][M]{};

  void init(const State &state) {
    memcpy(rho, state.rho, sizeof(rho));
    memcpy(p, state.p, sizeof(p));
    memcpy(old_p, state.old_p, sizeof(old_p));
    g = state.g;
    memcpy(field, state.field, sizeof(field));
    memcpy(dirs, state.dirs, sizeof(dirs));
    memcpy(last_use, state.last_use, sizeof(last_use));
    memcpy(velocity.v, state.velocity.v, sizeof(velocity.v));
    memcpy(velocity_flow.v, state.velocity_flow.v, sizeof(velocity_flow.v));
    UT = state.UT;
  }

  State state() const {
    auto ret = State{};
    memcpy(ret.rho, rho, sizeof(rho));
    memcpy(ret.p, p, sizeof(p));
    memcpy(ret.old_p, old_p, sizeof(old_p));
    ret.g = g;
    memcpy(ret.field, field, sizeof(field));
    memcpy(ret.dirs, dirs, sizeof(dirs));
    memcpy(ret.last_use, last_use, sizeof(last_use));
    memcpy(ret.velocity.v, velocity.v, sizeof(velocity.v));
    memcpy(ret.velocity_flow.v, velocity_flow.v, sizeof(velocity_flow.v));
    ret.UT = UT;
    return ret;
  }

  void step(size_t i) {
    PT total_delta_p = 0;
    // Apply external forces
    for (size_t x = 0; x < N; ++x) {
      for (size_t y = 0; y < M; ++y) {
        if (field[x][y] == '#')
          continue;
        if (field[x + 1][y] != '#')
          velocity.add(x, y, 1, 0, g);
      }
    }

    // Apply forces from p
    memcpy(old_p, p, sizeof(p));
    for (size_t x = 0; x < N; ++x) {
      for (size_t y = 0; y < M; ++y) {
        if (field[x][y] == '#')
          continue;
        for (auto [dx, dy] : deltas) {
          int nx = x + dx, ny = y + dy;
          if (field[nx][ny] != '#' && old_p[nx][ny] < old_p[x][y]) {
            auto delta_p = old_p[x][y] - old_p[nx][ny];
            auto force = delta_p;
            auto &contr = velocity.get(nx, ny, -dx, -dy);
            if (contr * to_vt(rho[(int)field[nx][ny]]) >= to_vt(force)) {
              contr -= to_vt(force / rho[(int)field[nx][ny]]);
              continue;
            }
            force -= to_pt(contr) * rho[(int)field[nx][ny]];
            contr = 0;
            velocity.add(x, y, dx, dy, to_vt(force / rho[(int)field[x][y]]));
            p[x][y] -= force / to_pt(dirs[x][y]);
            total_delta_p -= force / to_pt(dirs[x][y]);
          }
        }
      }
    }

    // Make flow from velocities
    velocity_flow = {};
    bool prop = false;
    do {
      UT += 2;
      prop = 0;
      for (size_t x = 0; x < N; ++x) {
        for (size_t y = 0; y < M; ++y) {
          if (field[x][y] != '#' && last_use[x][y] != UT) {
            auto [t, local_prop, _] = propagate_flow(x, y, 1);
            if (t > VT(0)) {
              prop = 1;
            }
          }
        }
      }
    } while (prop);

    // Recalculate p with kinetic energy
    for (size_t x = 0; x < N; ++x) {
      for (size_t y = 0; y < M; ++y) {
        if (field[x][y] == '#')
          continue;
        for (auto [dx, dy] : deltas) {
          auto old_v = velocity.get(x, y, dx, dy);
          auto new_v = velocity_flow.get(x, y, dx, dy);
          if (old_v > to_vt(0)) {
            assert(to_vt(new_v) <= old_v);
            velocity.get(x, y, dx, dy) = to_vt(new_v);
            auto force =
                (to_vft(old_v) - new_v) * to_vft(rho[(int)field[x][y]]);
            if (field[x][y] == '.')
              force *= to_vft(0.8);
            if (field[x + dx][y + dy] == '#') {
              p[x][y] += to_pt(force / to_vft(dirs[x][y]));
              total_delta_p += to_pt(force / to_vft(dirs[x][y]));
            } else {
              p[x + dx][y + dy] += to_pt(force / to_vft(dirs[x + dx][y + dy]));
              total_delta_p += to_pt(force / to_vft(dirs[x + dx][y + dy]));
            }
          }
        }
      }
    }

    UT += 2;
    prop = false;
    for (size_t x = 0; x < N; ++x) {
      for (size_t y = 0; y < M; ++y) {
        if (field[x][y] != '#' && last_use[x][y] != UT) {
          if (random01() < move_prob(x, y)) {
            prop = true;
            propagate_move(x, y, true);
          } else {
            propagate_stop(x, y, true);
          }
        }
      }
    }

    if (prop) {
      cout << "Tick " << i << ":\n";
      for (size_t x = 0; x < N; ++x) {
        cout << field[x] << "\n";
      }
    }
  }
};
