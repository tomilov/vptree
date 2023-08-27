#include <cmath>
#include <cstddef>

namespace
{

using F = double;

constexpr F kEps = F(1) / F(size_t(1) << 24);

struct Point
{
    F x, y;
    size_t i;

    F rho() const
    {
        return std::sqrt(rho2());
    }

    F rho2() const
    {
        return x * x + y * y;
    }

    Point& operator -= (const Point& rhs)
    {
        x -= rhs.x;
        y -= rhs.y;
        return *this;
    }

    F rho(Point point) const
    {
        point -= *this;
        return point.rho();
    }

    F rho2(Point point) const
    {
        point -= *this;
        return point.rho2();
    }

    bool closer(Point lhs, Point rhs) const
    {
        lhs -= *this;
        rhs -= *this;
        auto l = lhs.rho();
        auto r = rhs.rho();
        if (l + kEps < r) {
            return true;
        } else if (r + kEps < l) {
            return false;
        } else if (lhs.x + kEps < rhs.x) {
            return true;
        } else if (rhs.x + kEps < lhs.x) {
            return false;
        } else if (lhs.y + kEps < rhs.y) {
            return true;
        } else if (rhs.y + kEps < lhs.y) {
            return false;
        } else {
            return lhs.y * rhs.x < lhs.x * rhs.y;  // CCW
        }
    }
};

}

#include <fmt/color.h>
#include <fmt/format.h>
#include <fmt/chrono.h>

template<>
struct fmt::formatter<Point> : fmt::formatter<fmt::string_view>
{
    template<typename FormatContext>
    auto format(const Point& p, FormatContext& ctx) const
    {
        return fmt::format_to(ctx.out(), "Point{{.x = {:.3f}, .y = {:.3f}, .i = {}}}", p.x, p.y, p.i);
    }
};

#include <chrono>
#include <utility>

namespace
{

auto makeCheckpoint()
{
    return [start = std::chrono::high_resolution_clock::now()](fmt::string_view message = {}) mutable
    {
        const auto now = std::chrono::high_resolution_clock::now();
        const auto dt = std::chrono::ceil<std::chrono::milliseconds>(now - std::exchange(start, now));
        if (message.size() == 0) {
            fmt::println("{}", dt);
        } else {
            fmt::println("{}: {}", message, dt);
        }
        return dt;
    };
}

}

#include <algorithm>
#include <iterator>
#include <vector>
#include <optional>

#include <cassert>

namespace
{

template<typename I, typename IterLess>
std::pair<I, I> partition(I l, I mid, I r, const IterLess& iterLess)
{
    assert(l != r);
    std::iter_swap(mid, std::prev(r));
    mid = std::prev(r);
    for (auto i = l; i != mid; ++i) {
        if (iterLess(i, std::prev(r))) {
            std::iter_swap(i, l++);
        } else if (!iterLess(std::prev(r), i)) {
            if (--mid == i) {
                break;
            }
            std::iter_swap(mid, i);
        }
    }
    if (std::distance(l, mid) < std::distance(mid, r)) {
        return {l, std::rotate(l, mid, r)};
    } else {
        return {l, std::swap_ranges(mid, r, l)};
    }
}

template<typename I, typename IterLess>
void quickSelect(I l, const I mid, I r, const IterLess& iterLess)
{
    for (;;) {
        auto [ll, rr] = partition(l, mid, r, iterLess);
        if (std::distance(l, mid) < std::distance(l, ll)) {
            r = ll;
        } else if (std::distance(rr, r) < std::distance(mid, r)) {
            break;
        } else {
            l = rr;
        }
    }
}

template<typename I, typename C>
struct PointLocation
{
    explicit PointLocation(ptrdiff_t leafSize)
        : leafSize{leafSize}
    {
        assert(0 < leafSize);
    }

    // sites will be reordered
    std::vector<std::vector<C>> operator()(I sl, I sr, C ql, C qr)
    {
        assert(std::distance(sl, sr) > 1);
        assert(std::distance(ql, qr) > 0);
        assert(std::empty(vpTree));

        std::chrono::milliseconds total = {};

        {
            auto checkpoint = makeCheckpoint();
            auto vantagePoint = ql;
            if (!build(sl, sr, ql, qr, vantagePoint)) {
                vpTree.push_back({.vantagePoint = vantagePoint});
            }
            total += checkpoint(fmt::format(fmt::fg(fmt::color::blue), "build"));
        }

        fmt::println("vp tree size: {}", std::size(vpTree));

        std::vector<std::vector<C>> results;
        results.reserve(size_t(std::distance(ql, qr)));
        {
            auto checkpoint = makeCheckpoint();
            for (auto queryPoint = ql; queryPoint != qr; ++queryPoint) {
                searchClosest(sl, sr, *queryPoint, vpTree.back(), results.emplace_back());
            }
            total += checkpoint(fmt::format(fmt::fg(fmt::color::green), "query"));
        }

        fmt::println("total time: {}\n", fmt::styled(total, fmt::fg(fmt::color::red)));

        vpTree.clear();
        return results;
    }

private:
    struct Node
    {
        C vantagePoint;
        std::optional<size_t> l, r;
    };

    const ptrdiff_t leafSize;
    std::vector<Node> vpTree;

    std::optional<size_t> build(I sl, I sr, C vl, C vr, C& vantagePoint)
    {
        if (std::distance(sl, sr) <= leafSize) {
            return std::nullopt;
        }
        auto mid = std::next(sl, std::distance(sl, sr) / 2);
        quickSelect(sl, mid, sr, [&vantagePoint](C lhs, C rhs) { return vantagePoint->closer(*lhs, *rhs); });

        Node node = {
            .vantagePoint = vantagePoint,
        };
        if (++vantagePoint == vr) {
            vantagePoint = vl;
        }
        node.l = build(sl, mid, vl, vr, vantagePoint);
        node.r = build(std::next(mid), sr, vl, vr, vantagePoint);
        vpTree.push_back(node);
        return std::size(vpTree) - 1;
    }

    void tryAddResult(C s, const Point& queryPoint, std::vector<C>& result) const
    {
        if (std::empty(result)) {
            result.push_back(s);
        } else {
            auto lhs = queryPoint.rho(*result.front());
            auto rhs = queryPoint.rho(*s);
            if (!(lhs + kEps < rhs)) {
                if (rhs + kEps < lhs) {
                    result.clear();
                }
                result.push_back(s);
                if (std::size(result) > 1) {
                    const auto closer = [&queryPoint](C lhs, C rhs) -> bool
                    {
                        return queryPoint.rho2(*lhs) < queryPoint.rho2(*rhs);
                    };
                    std::push_heap(std::begin(result), std::end(result), closer);
                }
            }
        }
    }

    void searchClosestL(C sl, C sr, const Point& queryPoint, const Node& node, std::vector<C>& result) const
    {
        if (node.l) {
            return searchClosest(sl, sr, queryPoint, vpTree[*node.l], result);
        }
        for (auto s = sl; s != sr; ++s) {
            tryAddResult(s, queryPoint, result);
        }
    }

    void searchClosestR(C sl, C sr, const Point& queryPoint, const Node& node, std::vector<C>& result) const
    {
        if (node.r) {
            return searchClosest(sl, sr, queryPoint, vpTree[*node.r], result);
        }
        for (auto s = sl; s != sr; ++s) {
            tryAddResult(s, queryPoint, result);
        }
    }

    void searchClosest(C sl, C sr, const Point& queryPoint, const Node& node, std::vector<C>& result) const
    {
        auto mid = std::next(sl, std::distance(sl, sr) / 2);
        tryAddResult(mid, queryPoint, result);
        auto rho = node.vantagePoint->rho(queryPoint);
        auto threshold = node.vantagePoint->rho(*mid);
        bool left = rho + kEps < threshold;
        if (left) {
            searchClosestL(sl, mid, queryPoint, node, result);
        } else {
            searchClosestR(std::next(mid), sr, queryPoint, node, result);
        }
        if (!(queryPoint.rho(*result.front()) + kEps < std::abs(threshold - rho))) {
            if (left) {
                searchClosestR(std::next(mid), sr, queryPoint, node, result);
            } else {
                searchClosestL(sl, mid, queryPoint, node, result);
            }
        }
    }
};

}

#include <functional>
#include <numeric>
#include <random>
#include <numbers>

#include <cstdio>

int main()
{
    std::vector<Point> sites(100000);
    std::vector<Point> queryPoints(10000);

    {
        constexpr F kSize = F(10000);
        std::mt19937 g;
        auto checkpoint = makeCheckpoint();
        if ((false)) {
            std::uniform_real_distribution<F> r{-kSize, kSize};
            size_t m = 0;
            for (auto& [x, y, i]: sites) {
                x = r(g);
                y = r(g);
                i = m++;
            }
            size_t n = 0;
            for (auto& [x, y, i]: queryPoints) {
                x = r(g);
                y = r(g);
                i = n++;
            }
        } else {
            {
                const auto dPhi = F(2) * std::numbers::pi_v<F> / F(std::size(sites));
                size_t m = 0;
                for (auto& [x, y, i]: sites) {
                    i = m++;
                    auto phi = F(i) * dPhi;
                    x = kSize * std::cos(phi);
                    y = kSize * std::sin(phi);
                }
            }
            {
                constexpr F kRaius = kSize / F(994);
                size_t n = 0;
                const auto dPhi = F(2) * std::numbers::pi_v<F> / F(std::size(queryPoints));
                for (auto& [x, y, i]: queryPoints) {
                    i = n++;
                    auto phi = F(i) * dPhi;
                    x = kRaius * std::cos(phi);
                    y = kRaius * std::sin(phi);
                }
            }

            std::shuffle(std::begin(sites), std::end(sites), g);
            std::shuffle(std::begin(queryPoints), std::end(queryPoints), g);
        }
        checkpoint("gen input");
    }

    //fmt::println("{}", fmt::join(sites, ", "));
    //fmt::println("{}", fmt::join(queryPoints, ", "));

#if 1
    for (ptrdiff_t leafSize: {25}) {
#else
    for (ptrdiff_t leafSize : {1, 2, 4, 8, 16, 32, 64, 128, 256}) {
#endif
        fmt::println("leaf size: {}", leafSize);

        using I = typename std::vector<Point>::iterator;
        using C = typename std::vector<Point>::const_iterator;
        PointLocation<I, C> pointLocation{leafSize};
        auto results = pointLocation(std::begin(sites), std::end(sites), std::cbegin(queryPoints), std::cend(queryPoints));

        std::fflush(stdout);
        continue;

        auto orderIsLess = [](C lhs, C rhs) { return lhs->i < rhs->i; };
        auto q = std::cbegin(queryPoints);
        for (auto& result: results) {
            auto& queryPoint = *q++;
            const auto checkIsClose = [rho = queryPoint.rho(*result.front()), &queryPoint](C s) -> bool
            {
                return !(rho + kEps < queryPoint.rho(*s));
            };

            result.erase(std::partition(std::next(std::begin(result)), std::end(result), checkIsClose), std::end(result));
            std::sort(std::begin(result), std::end(result), orderIsLess);

            std::vector<C> expected;
            {
                expected.reserve(std::size(result));

                auto s = std::cbegin(sites);
                expected.push_back(s);
                auto lhs = queryPoint.rho(*s);
                while (++s != std::cend(sites)) {
                    auto rhs = queryPoint.rho(*s);
                    if (!(lhs + kEps < rhs)) {
                        if (rhs + kEps < lhs) {
                            expected.clear();
                            lhs = rhs;
                        }
                        expected.push_back(s);
                        if (std::size(expected) > 1) {
                            const auto closer = [&queryPoint](C lhs, C rhs) -> bool
                            {
                                return queryPoint.rho2(*lhs) < queryPoint.rho2(*rhs);
                            };
                            std::push_heap(std::begin(expected), std::end(expected), closer);
                            lhs = queryPoint.rho(*expected.front());
                        }
                    }
                }
            }
            result.erase(std::partition(std::next(std::begin(expected)), std::end(expected), checkIsClose), std::end(expected));
            std::sort(std::begin(expected), std::end(expected), orderIsLess);

            if (!std::equal(std::cbegin(result), std::cend(result), std::cbegin(expected), std::cend(expected), [](C lhs, C rhs) { return lhs->i == rhs->i; })) {
                fmt::println("Query point\n\t{}", queryPoint);
                for (auto site: result) {
                    fmt::println("\tr {}: {}", *site, queryPoint.rho(*site));
                }
                for (const auto& e: expected) {
                    fmt::println("\te {}: {}", *e, queryPoint.rho(*e));
                }
                fmt::println("");
            }
        }
    }
}
