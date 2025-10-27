<?xml version='1.0' encoding='UTF-8' standalone='yes' ?>
<tagfile doxygen_version="1.12.0">
  <compound kind="file">
    <name>blocked.hpp</name>
    <path>mumosa/</path>
    <filename>blocked_8hpp.html</filename>
    <includes id="simple_8hpp" name="simple.hpp" local="yes" import="no" module="no" objc="no">simple.hpp</includes>
    <class kind="struct">mumosa::BlockedOptions</class>
    <class kind="struct">mumosa::BlockedWorkspace</class>
    <class kind="class">mumosa::BlockedIndicesFactory</class>
    <class kind="struct">mumosa::BlockedIndicesFactory::Buffers</class>
    <namespace>mumosa</namespace>
  </compound>
  <compound kind="file">
    <name>combine_scaled_embeddings.hpp</name>
    <path>mumosa/</path>
    <filename>combine__scaled__embeddings_8hpp.html</filename>
    <namespace>mumosa</namespace>
  </compound>
  <compound kind="file">
    <name>compute_scale.hpp</name>
    <path>mumosa/</path>
    <filename>compute__scale_8hpp.html</filename>
    <namespace>mumosa</namespace>
  </compound>
  <compound kind="file">
    <name>mumosa.hpp</name>
    <path>mumosa/</path>
    <filename>mumosa_8hpp.html</filename>
    <includes id="simple_8hpp" name="simple.hpp" local="yes" import="no" module="no" objc="no">simple.hpp</includes>
    <includes id="blocked_8hpp" name="blocked.hpp" local="yes" import="no" module="no" objc="no">blocked.hpp</includes>
    <includes id="compute__scale_8hpp" name="compute_scale.hpp" local="yes" import="no" module="no" objc="no">compute_scale.hpp</includes>
    <includes id="combine__scaled__embeddings_8hpp" name="combine_scaled_embeddings.hpp" local="yes" import="no" module="no" objc="no">combine_scaled_embeddings.hpp</includes>
    <namespace>mumosa</namespace>
  </compound>
  <compound kind="file">
    <name>simple.hpp</name>
    <path>mumosa/</path>
    <filename>simple_8hpp.html</filename>
    <class kind="struct">mumosa::Options</class>
    <namespace>mumosa</namespace>
  </compound>
  <compound kind="class">
    <name>mumosa::BlockedIndicesFactory</name>
    <filename>classmumosa_1_1BlockedIndicesFactory.html</filename>
    <templarg>typename Index_</templarg>
    <templarg>typename Block_</templarg>
    <class kind="struct">mumosa::BlockedIndicesFactory::Buffers</class>
    <member kind="function">
      <type></type>
      <name>BlockedIndicesFactory</name>
      <anchorfile>classmumosa_1_1BlockedIndicesFactory.html</anchorfile>
      <anchor>a0443db4d10ca9b69e2f1999be3d9a743</anchor>
      <arglist>(const Index_ num_cells, const Block_ *block)</arglist>
    </member>
    <member kind="function">
      <type>const std::vector&lt; Index_ &gt; &amp;</type>
      <name>sizes</name>
      <anchorfile>classmumosa_1_1BlockedIndicesFactory.html</anchorfile>
      <anchor>a8825d775eb735f72c5c0be03ed860a66</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>Buffers&lt; Input_ &gt;</type>
      <name>create_buffers</name>
      <anchorfile>classmumosa_1_1BlockedIndicesFactory.html</anchorfile>
      <anchor>a94a2caff1b1b77b7e6ca370cb1d57970</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>build</name>
      <anchorfile>classmumosa_1_1BlockedIndicesFactory.html</anchorfile>
      <anchor>a14a3376745471edd1b72e6c508d56b49</anchor>
      <arglist>(const std::size_t num_dim, const Input_ *const data, const knncolle::Builder&lt; Index_, Input_, Distance_, Matrix_ &gt; &amp;builder, std::vector&lt; std::shared_ptr&lt; const knncolle::Prebuilt&lt; Index_, Input_, Distance_ &gt; &gt; &gt; &amp;output, Buffers&lt; Input_ &gt; &amp;work) const</arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; std::shared_ptr&lt; const knncolle::Prebuilt&lt; Index_, Input_, Distance_ &gt; &gt; &gt;</type>
      <name>build</name>
      <anchorfile>classmumosa_1_1BlockedIndicesFactory.html</anchorfile>
      <anchor>ab7bfb8020b5379643ab1bd9b0fab806e</anchor>
      <arglist>(const std::size_t num_dim, const Input_ *const data, const knncolle::Builder&lt; Index_, Input_, Distance_, Matrix_ &gt; &amp;builder) const</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>mumosa::BlockedOptions</name>
    <filename>structmumosa_1_1BlockedOptions.html</filename>
    <member kind="variable">
      <type>int</type>
      <name>num_neighbors</name>
      <anchorfile>structmumosa_1_1BlockedOptions.html</anchorfile>
      <anchor>af5844d64ee5524d9074f97a8d6cbc237</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>scran_blocks::WeightPolicy</type>
      <name>block_weight_policy</name>
      <anchorfile>structmumosa_1_1BlockedOptions.html</anchorfile>
      <anchor>aa39215f15f49244c734e5d6fdde99559</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>scran_blocks::VariableWeightParameters</type>
      <name>variable_block_weight_parameters</name>
      <anchorfile>structmumosa_1_1BlockedOptions.html</anchorfile>
      <anchor>a3da1b340177201902deb499edf632218</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>num_threads</name>
      <anchorfile>structmumosa_1_1BlockedOptions.html</anchorfile>
      <anchor>abf104360e28bdb625b49ab062d9066bc</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>mumosa::BlockedWorkspace</name>
    <filename>structmumosa_1_1BlockedWorkspace.html</filename>
    <templarg>typename Distance_</templarg>
  </compound>
  <compound kind="struct">
    <name>mumosa::BlockedIndicesFactory::Buffers</name>
    <filename>structmumosa_1_1BlockedIndicesFactory_1_1Buffers.html</filename>
    <templarg>typename Input_</templarg>
  </compound>
  <compound kind="struct">
    <name>mumosa::Options</name>
    <filename>structmumosa_1_1Options.html</filename>
    <member kind="variable">
      <type>int</type>
      <name>num_neighbors</name>
      <anchorfile>structmumosa_1_1Options.html</anchorfile>
      <anchor>a76621d12b7ed8b1c072ccf3598761640</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>num_threads</name>
      <anchorfile>structmumosa_1_1Options.html</anchorfile>
      <anchor>a390a953ace5a6fa0a172efb8cb691929</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="namespace">
    <name>mumosa</name>
    <filename>namespacemumosa.html</filename>
    <class kind="class">mumosa::BlockedIndicesFactory</class>
    <class kind="struct">mumosa::BlockedOptions</class>
    <class kind="struct">mumosa::BlockedWorkspace</class>
    <class kind="struct">mumosa::Options</class>
    <member kind="function">
      <type>BlockedWorkspace&lt; Distance_ &gt;</type>
      <name>create_workspace</name>
      <anchorfile>namespacemumosa.html</anchorfile>
      <anchor>ad4d80dfac74448a389f91f72f8a31548</anchor>
      <arglist>(const std::vector&lt; Index_ &gt; &amp;block_sizes, const BlockedOptions &amp;options)</arglist>
    </member>
    <member kind="function">
      <type>std::pair&lt; Distance_, Distance_ &gt;</type>
      <name>compute_distance_blocked</name>
      <anchorfile>namespacemumosa.html</anchorfile>
      <anchor>a865e2bcb5ff9ae07f7f4d9c762720f29</anchor>
      <arglist>(const std::vector&lt; std::shared_ptr&lt; const knncolle::Prebuilt&lt; Index_, Input_, Distance_ &gt; &gt; &gt; &amp;prebuilts, BlockedWorkspace&lt; Distance_ &gt; &amp;workspace, const BlockedOptions &amp;options)</arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; std::shared_ptr&lt; const knncolle::Prebuilt&lt; Index_, Input_, Distance_ &gt; &gt; &gt;</type>
      <name>build_blocked_indices</name>
      <anchorfile>namespacemumosa.html</anchorfile>
      <anchor>abd487dda75494406f5e8b66e5d2d68f1</anchor>
      <arglist>(const std::size_t num_dim, const std::vector&lt; Index_ &gt; block_sizes, const Input_ *const data, const knncolle::Builder&lt; Index_, Input_, Distance_, Matrix_ &gt; &amp;builder)</arglist>
    </member>
    <member kind="function">
      <type>std::pair&lt; Distance_, Distance_ &gt;</type>
      <name>compute_distance_blocked</name>
      <anchorfile>namespacemumosa.html</anchorfile>
      <anchor>a91388ad0732cac79997e01d9b9875372</anchor>
      <arglist>(const std::size_t num_dim, const std::vector&lt; Index_ &gt; &amp;block_sizes, const Input_ *const data, const knncolle::Builder&lt; Index_, Input_, Distance_, Matrix_ &gt; &amp;builder, const BlockedOptions &amp;options)</arglist>
    </member>
    <member kind="function">
      <type>std::pair&lt; Distance_, Distance_ &gt;</type>
      <name>compute_distance_blocked</name>
      <anchorfile>namespacemumosa.html</anchorfile>
      <anchor>a8510a19aca6de4983bcfbef0b97cf940</anchor>
      <arglist>(const std::size_t num_dim, const Index_ num_cells, const Input_ *const data, const Block_ *const block, const knncolle::Builder&lt; Index_, Input_, Distance_, Matrix_ &gt; &amp;builder, const BlockedOptions &amp;options)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>combine_scaled_embeddings</name>
      <anchorfile>namespacemumosa.html</anchorfile>
      <anchor>adc1aa9b1eba884ed114f58fdb9fd3b19</anchor>
      <arglist>(const std::vector&lt; std::size_t &gt; &amp;num_dims, const Index_ num_cells, const std::vector&lt; Input_ * &gt; &amp;embeddings, const std::vector&lt; Scale_ &gt; &amp;scaling, Output_ *const output)</arglist>
    </member>
    <member kind="function">
      <type>Distance_</type>
      <name>compute_scale</name>
      <anchorfile>namespacemumosa.html</anchorfile>
      <anchor>a8227367fe83f0bec0c50156c5f24a3b7</anchor>
      <arglist>(const std::pair&lt; Distance_, Distance_ &gt; &amp;ref, const std::pair&lt; Distance_, Distance_ &gt; &amp;target)</arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; Distance_ &gt;</type>
      <name>compute_scale</name>
      <anchorfile>namespacemumosa.html</anchorfile>
      <anchor>a26c6d69157bee4c3c7d595570761793b</anchor>
      <arglist>(const std::vector&lt; std::pair&lt; Distance_, Distance_ &gt; &gt; &amp;distances)</arglist>
    </member>
    <member kind="function">
      <type>std::pair&lt; Distance_, Distance_ &gt;</type>
      <name>compute_distance</name>
      <anchorfile>namespacemumosa.html</anchorfile>
      <anchor>ab9bb10478b865f21008d5b87e5d8cf4a</anchor>
      <arglist>(const Index_ num_cells, Distance_ *const distances)</arglist>
    </member>
    <member kind="function">
      <type>std::pair&lt; Distance_, Distance_ &gt;</type>
      <name>compute_distance</name>
      <anchorfile>namespacemumosa.html</anchorfile>
      <anchor>aae9909a1adc9470d2689046b7dab20df</anchor>
      <arglist>(const knncolle::Prebuilt&lt; Index_, Input_, Distance_ &gt; &amp;prebuilt, Distance_ *const distances, const Options &amp;options)</arglist>
    </member>
    <member kind="function">
      <type>std::pair&lt; Distance_, Distance_ &gt;</type>
      <name>compute_distance</name>
      <anchorfile>namespacemumosa.html</anchorfile>
      <anchor>ad5544f59301aea3947788e64190009d1</anchor>
      <arglist>(const std::size_t num_dim, const Index_ num_cells, const Input_ *const data, const knncolle::Builder&lt; Index_, Input_, Distance_, Matrix_ &gt; &amp;builder, const Options &amp;options)</arglist>
    </member>
  </compound>
  <compound kind="page">
    <name>index</name>
    <title>Multi-modal single-cell analyses</title>
    <filename>index.html</filename>
    <docanchor file="index.html">md__2github_2workspace_2README</docanchor>
  </compound>
</tagfile>
