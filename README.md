# PeakSegOptimal - Unconstrained FPOP Implementation

## Summary -

Implemented an unconstrained version of the FPOP algorithm for optimal changepoint detection with Poisson loss, extending the existing PeakSegOptimal package.

## Implementation

### Files Modified/Created:

1. **`src/funPieceListLog.h`** - Added `set_to_unconstrained_min_of` declaration
2. **`src/funPieceListLog.cpp`** - Implemented `set_to_unconstrained_min_of` method
3. **`src/PeakSegFPOPLogUnconstrained.cpp`** - Main algorithm implementation
4. **`src/interface.cpp`** - Added C interface wrapper
5. **`R/PeakSegFPOPLogUnconstrained.R`** - R wrapper function
6. **`tests/testthat/test-unconstrained.R`** - Test suite

### Key Changes:

- Implemented `PiecewisePoissonLossLog::set_to_unconstrained_min_of()` which computes the minimum cost function without up/down constraints
- Created `PeakSegFPOPLogUnconstrained()` algorithm that uses N×1 cost matrix (vs N×2 for constrained version)
- Added comprehensive tests comparing with `Segmentor3IsBack::Segmentor(model=1)`

## Testing

All tests pass, verifying that the unconstrained implementation:
- Correctly detects changepoints for varying penalty parameters
- Matches Segmentor3IsBack results
- Handles edge cases (constant data, distinct groups)
- Shows proper penalty-segment count relationship

## Links

- **Pull Request:** [link when you create it]
- **Test Results:** [link to test output]