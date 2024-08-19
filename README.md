# Multi-modal single-cell analyses

![Unit tests](https://github.com/libscran/mumosa/actions/workflows/run-tests.yaml/badge.svg)
![Documentation](https://github.com/libscran/mumosa/actions/workflows/doxygenate.yaml/badge.svg)
[![Codecov](https://codecov.io/gh/libscran/mumosa/graph/badge.svg?token=6lt7ykKjoH)](https://codecov.io/gh/libscran/mumosa)

## Overview

In multi-modal single-cell experiments, we obtain data of different modalities (e.g., RNA, protein) from the same set of cells.
Naturally, we would like to combine data from different modalities to increase the information available for each cell.
This is most relevant to analysis steps that yield results in terms of cells, like clustering and visualization with t-SNE or UMAP.
A simple combining strategy is to just concatenate the per-modality data matrices together into a single matrix for further analysis.
While convenient and compatible with many downstream procedures, this is complicated by the differences in the variance between modalities. 
Higher noise in one modality might drown out biological signal in another modality that has lower variance.

The **mumosa** algorithm scales embeddings to equalize noise across modalities prior to concatenation.
First, we compute the median distance to the $k$-th nearest neighbor in the low-dimensional embedding for each modality (e.g., after PCA).
This distance is used as a proxy for the modality-specific noise within a subpopulation containing at least $k$ cells.
We define a scaling factor for each modality as the ratio of the median distances for that modality compared to a "reference" modality.
We scale the modality's embedding by its factor, removing differences in variance due to irrelevant factors like the scale of expression values, number of features, etc.
We then concatenate the matrices to form a single embedding for further analysis.

## Quick start

Given the embedding coordinates for multiple modalities, we compute the median distance to the $k$-nearest neighbor for each modality: 

```cpp
#include "mumosa/mumosa.hpp"

// Mocking up some modalities.
int nobs = 1000;
std::vector<int> dimensions(3, 20);
std::vector<std::vector<double> > embeddings(3);
for (int m = 0; m < 3; ++m) {
    embeddings[m].resize(nobs * dimensions[m]);
}

// Computing distances per modality.
mumosa::Options opt;
opt.num_neighbors = 20;
opt.num_threads = 3;

std::vector<std::pair<double, double> > distances;
for (int m = 0; m < 3; ++m) {
    distances[m] = mumosa::compute_distance(
        dimensions[m],
        nobs,
        embeddings[m].data(),
        knncolle::VptreeBuilder<>(), // any NN algorithm can be used here.
        opt
    );
}
```

We then compute scaling factors for each modality:

```cpp
auto scale = mumosa::compute_scale(distances);
```

And combine the scaled per-modality embeddings into a single matrix:

```cpp
size_t ntotal = std::accumulate(dimensions.begin(), dimensions.end(), 0);
std::vector<double> combined(ntotal * nobs);

std::vector<const double*> inputs;
for (const auto& em : embeddings) {
    inputs.push_back(em.data());
}

mumosa::combine_scaled_embeddings(
    dimensions,
    nobs,
    inputs,
    scale,
    combined.data()
);
```

Check out the [reference documentation](https://libscran.github.io/mumosa) for more details.

## Further comments

The key element of the **mumosa** approach is the use of the distance to the nearest neighbors as a measure of (uninteresting) spread.
We do not use the total variance for each embedding as this includes the biological heterogeneity of interest.
Scaling by the total variance would reduce the contribution of the most relevant modalities, which is obviously not desirable.
**mumosa** aims to remove differences in the magnitude of noise while preserving modality-specific biological signal in the concatenated matrix.

The other appealing aspect of **mumosa** lies in its simplicity relative to other approaches (e.g., multi-modal factor analyses, intersection of simplicial sets). 
It returns a combined matrix of embeddings that can be directly used in any downstream analysis steps like clustering, t-SNE, UMAP, etc. without any issues.
No further transformations beyond scaling are performed, ensuring that population structure within each modality is faithfully represented in the combined embedding.
Most importantly, **mumosa** is very easy to implement, and that's probably what I like the most about it - keep it simple, stupid. 

Alright, now for some of the caveats.
The most obvious one is that the distance to a neighbor may not be an accurate relative measure of the within-population variance.
Even in the simplest cases of i.i.d. noise, the distance is not proportional to the standard deviation at lower dimensions 
(see analysis [here](tests/R/dimensions.Rmd)), to say nothing of heteroskedasticity across dimensions. 
Nonetheless, **mumosa** can still be useful for downstream procedures that perform distance calculations between cells,
as it ensures that each modality contributes equally to the distance between cells from the same subpopulation in the combined embedding.

A more subtle issue with the **mumosa** approach is that it assumes that the density of cells in each subpopulation is comparable across modalities.
If one modality has all of its cells in a single subpopulation, the median distance to a nearest neighbor will be lower, causing the modality's embedding to be scaled up unnecessarily.
This is a fundamentally difficult issue to fix as it requires a decision on what a "subpopulation" actually is.
For example, a big blob of cells may contain further interesting structure, in which case **mumosa**'s upscaling would be appropriate.
Users who know better (e.g., from control data) can adjust the scaling factors to give appropriate weights to each modality.

## Building projects

### CMake with `FetchContent`

If you're using CMake, you just need to add something like this to your `CMakeLists.txt`:

```cmake
include(FetchContent)

FetchContent_Declare(
  mumosa
  GIT_REPOSITORY https://github.com/libscran/mumosa
  GIT_TAG master # or any version of interest
)

FetchContent_MakeAvailable(mumosa)
```

Then you can link to **mumosa** to make the headers available during compilation:

```cmake
# For executables:
target_link_libraries(myexe libscran::mumosa)

# For libaries
target_link_libraries(mylib INTERFACE libscran::mumosa)
```

### CMake with `find_package()`

```cmake
find_package(libscran_mumosa CONFIG REQUIRED)
target_link_libraries(mylib INTERFACE libscran::mumosa)
```

To install the library, use:

```sh
mkdir build && cd build
cmake .. -DMUMOSA_TESTS=OFF
cmake --build . --target install
```

By default, this will use `FetchContent` to fetch all external dependencies.
If you want to install them manually, use `-DMUMOSA_FETCH_EXTERN=OFF`.
See the tags in [`extern/CMakeLists.txt`](extern/CMakeLists.txt) to find compatible versions of each dependency.

### Manual

If you're not using CMake, the simple approach is to just copy the files in `include/` - either directly or with Git submodules - and include their path during compilation with, e.g., GCC's `-I`.
This requires the external dependencies listed in [`extern/CMakeLists.txt`](extern/CMakeLists.txt), which also need to be made available during compilation.
