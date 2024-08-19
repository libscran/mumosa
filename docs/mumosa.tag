<?xml version='1.0' encoding='UTF-8' standalone='yes' ?>
<tagfile doxygen_version="1.9.8">
  <compound kind="file">
    <name>mumosa.hpp</name>
    <path>mumosa/</path>
    <filename>mumosa_8hpp.html</filename>
    <class kind="struct">mumosa::Options</class>
    <namespace>mumosa</namespace>
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
    <class kind="struct">mumosa::Options</class>
    <member kind="function">
      <type>std::pair&lt; Float_, Float_ &gt;</type>
      <name>compute_distance</name>
      <anchorfile>namespacemumosa.html</anchorfile>
      <anchor>a1dd47f25008fcc2d8c1249761a03b466</anchor>
      <arglist>(Index_ num_cells, Float_ *distances)</arglist>
    </member>
    <member kind="function">
      <type>std::pair&lt; Float_, Float_ &gt;</type>
      <name>compute_distance</name>
      <anchorfile>namespacemumosa.html</anchorfile>
      <anchor>a8bdd5bb9377ce3e6852eaae388921b89</anchor>
      <arglist>(const knncolle::Prebuilt&lt; Dim_, Index_, Float_ &gt; &amp;prebuilt, const Options &amp;options)</arglist>
    </member>
    <member kind="function">
      <type>std::pair&lt; Float_, Float_ &gt;</type>
      <name>compute_distance</name>
      <anchorfile>namespacemumosa.html</anchorfile>
      <anchor>a1759ae6e54125cf6286bb4db2f17eae4</anchor>
      <arglist>(Dim_ num_dim, Index_ num_cells, const Float_ *data, const knncolle::Builder&lt; knncolle::SimpleMatrix&lt; Dim_, Index_, Float_ &gt;, Float_ &gt; &amp;builder, const Options &amp;options)</arglist>
    </member>
    <member kind="function">
      <type>Float_</type>
      <name>compute_scale</name>
      <anchorfile>namespacemumosa.html</anchorfile>
      <anchor>a889cc7f4d99446e2a7dcd1d398763ace</anchor>
      <arglist>(const std::pair&lt; Float_, Float_ &gt; &amp;ref, const std::pair&lt; Float_, Float_ &gt; &amp;target)</arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; Float_ &gt;</type>
      <name>compute_scale</name>
      <anchorfile>namespacemumosa.html</anchorfile>
      <anchor>a9640853d15e0d836bfdf4d54a0874ac5</anchor>
      <arglist>(const std::vector&lt; std::pair&lt; Float_, Float_ &gt; &gt; &amp;distances)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>combine_scaled_embeddings</name>
      <anchorfile>namespacemumosa.html</anchorfile>
      <anchor>a5abe395d9107cd90aa10bba4362f847c</anchor>
      <arglist>(const std::vector&lt; Dim_ &gt; &amp;num_dims, Index_ num_cells, const std::vector&lt; Input_ * &gt; &amp;embeddings, const std::vector&lt; Scale_ &gt; &amp;scaling, Output_ *output)</arglist>
    </member>
  </compound>
  <compound kind="page">
    <name>index</name>
    <title>Multi-modal single-cell analyses</title>
    <filename>index.html</filename>
    <docanchor file="index.html">md__2github_2workspace_2README</docanchor>
  </compound>
</tagfile>
