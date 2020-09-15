model TwoElectrodeStimulus "Stimulus Created by two electrodes"
  import pi = Modelica.Constants.pi;
  import SI = Modelica.SIunits;
  parameter SI.Length x_anode = 100;
  parameter SI.Length y_anode = 100;
  parameter SI.Length x_cathode = 0.000;
  parameter SI.Length y_cathode = 1;
  input SI.Current I;
  parameter NeuronModelParameters params;
  output SI.Voltage V[params.Nnodes];
  SI.Length node_x[params.Nnodes];
  Real ra[params.Nnodes];
  Real rc[params.Nnodes];
  Real anode_potential[params.Nnodes];
  Real cathode_potential[params.Nnodes];
  Real epot[params.Nnodes];
  Real axon_length;
  Real axon_start;
equation
  axon_length = (params.Nnodes - 1) * params.internodal_distance;
  axon_start = -axon_length / 2;
  for i in 1:params.Nnodes loop
    node_x[i] = axon_start + (i - 1) * params.internodal_distance;
    ra[i] = sqrt((node_x[i] - x_anode) ^ 2 + y_anode ^ 2);
    rc[i] = sqrt((node_x[i] - x_cathode) ^ 2 + y_cathode ^ 2);
    anode_potential[i] = params.exoplasm_resistivity / (2 * pi * ra[i]);
    cathode_potential[i] = -params.exoplasm_resistivity / (2 * pi * rc[i]);
  end for;
  epot = anode_potential + cathode_potential;
  V = epot * I;
  annotation(Diagram(coordinateSystem(extent = {{-150, -90}, {150, 90}}, preserveAspectRatio = true, initialScale = 0.1, grid = {5, 5})), Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(visible = true, lineColor = {0, 114, 195}, fillColor = {255, 255, 255}, extent = {{-100, -100}, {100, 100}}, radius = 25), Text(visible = true, textColor = {64, 64, 64}, extent = {{-150, 110}, {150, 150}}, textString = "%name")}));
end TwoElectrodeStimulus;
