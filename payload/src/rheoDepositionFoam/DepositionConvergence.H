#ifndef DepositionConvergence_H
#define DepositionConvergence_H

#include <cmath>

// Pure numerical policy, also compiled by tests/test_convergence_policy.py.
namespace depositionConvergence
{
inline bool withinTolerance
(
    const double change, const double scale,
    const double relativeTolerance, const double absoluteTolerance
)
{
    return std::isfinite(change) && std::isfinite(scale)
        && change >= 0 && scale >= 0
        && change <= absoluteTolerance + relativeTolerance*scale;
}

class Tracker
{
    unsigned count_;
    double since_;
    double lastTime_;

public:
    Tracker() : count_(0), since_(-1), lastTime_(-1) {}
    void reset() { count_ = 0; since_ = -1; lastTime_ = -1; }
    unsigned count() const { return count_; }
    double duration(const double time) const
    {
        return since_ < 0 ? 0 : time-since_;
    }
    bool observe
    (
        const bool stable, const double time,
        const unsigned required, const double minimumTime,
        const double window
    )
    {
        // Repeated writes at one time are not independent convergence samples.
        if (!std::isfinite(time) || time <= lastTime_) return false;
        lastTime_ = time;
        if (!stable || time < minimumTime)
        {
            count_ = 0;
            since_ = -1;
            return false;
        }
        if (count_ == 0) since_ = time;
        ++count_;
        return count_ >= required && duration(time)+1e-12 >= window;
    }
};
}
#endif
