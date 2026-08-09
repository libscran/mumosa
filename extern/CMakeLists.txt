include(FetchContent)

FetchContent_Declare(
  knncolle
  GIT_REPOSITORY https://github.com/knncolle/knncolle
  GIT_TAG master # ^3.1.0
)

FetchContent_Declare(
  quickstats
  GIT_REPOSITORY https://github.com/LTLA/quickstats
  GIT_TAG master # ^0.2.0
)

FetchContent_Declare(
  scran_blocks
  GIT_REPOSITORY https://github.com/libscran/scran_blocks
  GIT_TAG master # ^0.1.0
)

FetchContent_Declare(
  sanisizer
  GIT_REPOSITORY https://github.com/LTLA/sanisizer
  GIT_TAG master # ^0.2.0
)

FetchContent_MakeAvailable(knncolle)
FetchContent_MakeAvailable(quickstats)
FetchContent_MakeAvailable(sanisizer)
FetchContent_MakeAvailable(scran_blocks)
