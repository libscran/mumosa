# Multi-modal single-cell analyses

![Unit tests](https://github.com/libscran/mumosa/actions/workflows/run-tests.yaml/badge.svg)
![Documentation](https://github.com/libscran/mumosa/actions/workflows/doxygenate.yaml/badge.svg)
[![Codecov](https://codecov.io/gh/libscran/mumosa/graph/badge.svg?token=6lt7ykKjoH)](https://codecov.io/gh/libscran/mumosa)

## Overview

Multi-modal single-cell experiments generate data of different modalities (e.g., RNA, protein) from the same set of cells.
In some parts of the analysis, we want to combine data from different modalities to increase the information available for each cell.
This is most relevant to results that are expressed in terms of cells, like clustering and visualization with t-SNE or UMAP.
A simple combining strategy is to just concatenate the per-modality data matrices together into a single matrix for further analysis.
While convenient and compatible with many downstream procedures, this is complicated by the differences in the variance between modalities - 
higher noise in one modality might drown out biological signal in another modality that has lower variance.

In **mumosa**, we scale embeddings to equalize noise across modalities prior to concatenation.
First, we compute the median distance to the $k$-th nearest neighbor in the low-dimensional embedding for each modality (e.g., after PCA).
This distance is used as a proxy for the modality-specific noise within a subpopulation containing at least $k$ cells.
We define a scaling factor for each modality as the ratio of the median distances for that modality compared to a "reference" modality.
We scale the modality's embedding by its factor, removing differences in variance due to irrelevant factors like the scale of expression values, number of features, etc.
We then concatenate the matrices to form a single embedding for further analysis.

We deliberately use the distance to the nearest neighbors instead of the total variance for each embedding.
The latter would unnecessarily penalize modalities with strong biological heterogeneity.
Our approach aims to remove differences in the magnitude of noise while preserving modality-specific biological signal in the concatenation.

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

We can then compute scaling factors for each modality:

```cpp
auto scale = mumosa::compute_scale(distances);
```

And use them to combine the embeddings into a single matrix:

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
cmake .. -DNENESUB_TESTS=OFF
cmake --build . --target install
```

By default, this will use `FetchContent` to fetch all external dependencies.
If you want to install them manually, use `-DNENESUB_FETCH_EXTERN=OFF`.
See the tags in [`extern/CMakeLists.txt`](extern/CMakeLists.txt) to find compatible versions of each dependency.

### Manual

If you're not using CMake, the simple approach is to just copy the files in `include/` - either directly or with Git submodules - and include their path during compilation with, e.g., GCC's `-I`.
This requires the external dependencies listed in [`extern/CMakeLists.txt`](extern/CMakeLists.txt), which also need to be made available during compilation.
