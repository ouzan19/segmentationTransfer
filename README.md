
<html>
<body>

<h1> Segmentation Transfer    </h1>

<h2> Compilation and usage:   </h2>

<p>  Use <b>CMakeLists</b>  to create Visual Studio 2013 solution. More commented and object-oriented code will be added soon.</p>

<p>  The user interface is done through key presses.</p>


<h2>  In <i>build\Release\faces</i> folder, there are 5 types of files:   </h2>


<ul style="list-style-type:disc">
  <li> .off files with "face" prefix are the inputs of the application </li>
  <li> .off files with "coarse" prefix are the output of PCA alignment method </li>
  <li> .off files with "fine" prefix are the output of ICP alignment method </li>
   <li> .off files with "deformed" prefix are the output of Laplacian alignment method </li>
    <li> .txt files with "corr" prefix are the corresponding indexes with reference model (face1.off) </li>
</ul>


<h2>Prerequisites</h2>

<ul style="list-style-type:disc">
  <li> VTK  <a href="http://www.vtk.org">  http://www.vtk.org</a> </li>
  <li> Eigen (included in source) </li>
</ul>

</body>
</html>