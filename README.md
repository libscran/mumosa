# Multi-modal single-cell analyses

![Unit tests](https://github.com/libscran/mumosa/actions/workflows/run-tests.yaml/badge.svg)
![Documentation](https://github.com/libscran/mumosa/actions/workflows/doxygenate.yaml/badge.svg)
[![Codecov](https://codecov.io/gh/libscran/mumosa/graph/badge.svg?token=6lt7ykKjoH)](https://codecov.io/gh/libscran/mumosa)

## Overview

In multi-modal single-cell experiments, we obtain data of different modalities (e.g., RNA, protein) from the same set of cells.
Naturally, we would like to combine data from different modalities to increase the information available for each cell in further analyses.
This is most relevant to analysis steps that operate on cells, e.g., clustering, visualization with t-SNE or UMAP.
The simplest combining strategy is to just concatenate the per-modality data matrices together into a single matrix for further analysis.
While convenient and compatible with many downstream procedures, this is complicated by the differences in the variance between modalities. 
Higher noise in one modality might drown out biological signal in another modality that has lower variance.

The **mumosa** algorithm scales the data from each modality to equalize "uninteresting" noise to concatenation.
First, we compute the median distance to the $k$-th nearest neighbor across all cells for each modality. 
This distance is used to as a measure of the modality-specific variance within each cell's local neighborhood.
We define a scaling factor for each modality as the ratio of the median distances for that modality compared to a "reference" modality.
We scale the modality's coordinates by this factor, removing differences in variance due to irrelevant factors like the scale of expression values, dimensionality, etc.
We then concatenate data across modalities into a single matrix for further analysis.

## Quick start

Each modality should be represented as a low-dimensional embedding (e.g., after PCA) for more efficient neighbor searches.
Given the embedding coordinates for multiple modalities, we compute the median distance to the $k$-nearest neighbor for each modality: 

```cpp
#include "mumosa/mumosa.hpp"

// Mocking up some modalities. For each modality 'm', we have a column-major
// array of 'embeddings[m]' of size 'dimensions[m] * nobs'. 
int nobs = 1000;
std::vector<int> dimensions(3, 20);
std::vector<std::vector<double> > embeddings(3);
for (int m = 0; m < 3; ++m) {
    embeddings[m].resize(nobs * dimensions[m]);
}

// Configuring the neighbor search algorithm; here, we'll be using an exact
// search based on VP trees with a Euclidean distance metric.
knncolle::VptreeBuilder<int, double, double> vp_builder(
    std::make_shared<knncolle::EuclideanDistance<double, double> >()
);

// Computing distances per modality.
mumosa::Options opt;
opt.num_neighbors = 20;
opt.num_threads = 3;

std::vector<std::pair<double, double> > distances;
std::vector<double> distbuffer(nobs);
for (int m = 0; m < 3; ++m) {
    distances[m] = mumosa::compute_distance(
        dimensions[m],
        nobs,
        embeddings[m].data(),
        vp_builder,
        distbuffer.data(),
        opt
    );
}
```

We compute scaling factors for each modality:

```cpp
auto scale = mumosa::compute_scale(distances);
```

And combine the scaled per-modality embeddings into a single matrix, which can be used for downstream steps like k-means clustering:

```cpp
std::size_t ntotal = std::accumulate(dimensions.begin(), dimensions.end(), 0);
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

The premise of the **mumosa** approach is that the median distance to the $k$-nearest neighbor is a suitable measure of (uninteresting) variation within each modality.
By quantifying the spread of cells in each local neighborhood, we capture the effects of dimensionality, scale, etc. without much contribution from biological variance. 
Scaling each modality by its median distance removes differences in the magnitude of noise while preserving modality-specific biological signal in the concatenated matrix.
In contrast, the total variance for each embedding includes the biological heterogeneity of interest.
Scaling each modality by its total variance would reduce the contribution of the most informative modalities, which is obviously not desirable.

Ideally, the median distance-to-neighbor would serve as a proxy for the average variance within subpopulations of at least $k + 1$ cells.
This provides an intuitive rationale for scaling each modality to equalize the within-population variance.
However, this interpretation has several caveats:

- Each modality may have a different subpopulation structure.
  A modality with a small number of large subpopulations will have a lower median distance-to-neighbor than a modality with a large number of small subpopulations,
  even if the variance within each subpopulation is the same - this would result in inappropriate upscaling of the former.
  In practice, this is not too problematic as the definition of a "subpopulation" is so vague that it's hard to say that our scaling is obviously wrong.
  For example, a big blob of cells may contain further interesting structure, in which case **mumosa**'s upscaling would be appropriate.
  Users who know better (e.g., from control data) can adjust the scaling factors to give appropriate weights to each modality.
- The median distance-to-neighbor is not an accurate relative measure of the variance at lower dimensions.
  Even in the simplest cases of i.i.d. noise, the distance is not proportional to the standard deviation at lower dimensions (see analysis [here](tests/R/dimensions.Rmd)). 
  Nonetheless, **mumosa** can still be useful for downstream procedures that perform distance calculations between cells,
  as it ensures that each modality contributes equally to the distance between cells from the same subpopulation in the combined embedding.

One appeal of **mumosa** is its simplicity relative to other approaches, e.g., multi-modal factor analyses, intersection of nearest-neighbor sets. 
No further transformations beyond scaling are performed, ensuring that population structure within each modality is faithfully represented in the combined embedding.
It is very easy to implement and the result is directly compatible with any downstream analysis step that can operate on an embedding matrix.
In fact, we only care about the median distance so we could save even more time by only performing the neighbor search for a subset of cells.

<!--- 
By comparison, if we did some probabilistic graph intersection, we'd be limited to algorithms that operate on weighted graphs.
And not all of these algorithms have the same interpretation of the weights, e.g., the t-SNE and UMAP probabilities are not the same.
So that would be a headache to recocile.
-->

## Building projects

### CMake with `FetchContent`

If you're using CMake, you just need to add something like this to your `CMakeLists.txt`:

```cmake
include(FetchContent)

FetchContent_Declare(
  mumosa
  GIT_REPOSITORY https://github.com/libscran/mumosa
  GIT_TAG master # replace with a pinned release
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

By default, this will use `FetchContent` to fetch all external dependencies.
Applications should consider pinning versions of all dependencies - see [`extern/CMakeLists.txt`](extern/CMakeLists.txt) for suggested versions.
If you want to install them manually, use `-DMUMOSA_FETCH_EXTERN=OFF`.

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

Again, this will use `FetchContent` to retrieve dependencies, see comments above.

### Manual

If you're not using CMake, the simple approach is to just copy the files in `include/` - either directly or with Git submodules - and include their path during compilation with, e.g., GCC's `-I`.
This also requires the external dependencies listed in [`extern/CMakeLists.txt`](extern/CMakeLists.txt). 
