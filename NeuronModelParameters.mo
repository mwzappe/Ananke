final record NeuronModelParameters "Parameters for the neuron model"
  import SI = Modelica.SIunits;
  import pi = Modelica.Constants.pi;
  /* Conversion Constants */
  /*
  import R = Modelica.Constants.R;
  import F = Modelica.Constants.F;
  constant SI.AmountOfSubstance mol = 1;
  constant SI.AmountOfSubstance mmol = 1e-3;
  constant SI.Volume L = SI.Conversions.from_litre(1);
  constant SI.Length cm = 0.01;
  constant SI.Length um = 1e-6;
  constant SI.Area cm2 = SI.Conversions.from_cm2(1);
  constant SI.Capacitance uF = 1e-6;
  constant SI.Resistance Ohm = 1;
  constant SI.Conductance mS = 1e-3;
  constant SI.Temperature T = 295.16;
  constant SI.Voltage mV = 1e-3;
*/
  // Reset units to SENN defaults
  constant Real R = 8.3144;
  constant Real F = 96.487;
  constant SI.AmountOfSubstance mol = 1e3;
  constant SI.AmountOfSubstance mmol = 1;
  constant SI.Volume L = 1;
  constant SI.Length cm = 1;
  constant SI.Length um = 1e-4;
  constant SI.Area cm2 = 1;
  constant SI.Capacitance uF = 1;
  constant SI.Resistance Ohm = 1;
  constant SI.Conductance mS = 1;
  constant SI.Temperature T = 295.16;
  constant SI.Voltage mV = 1;
  // Diameter of the nerve fiber
  parameter SI.Length fiber_diameter = 0.002 * cm;
  // Resitivity of the axoplasm
  parameter SI.Resistivity axoplasm_resistivity = 110.0 * Ohm * cm;
  // Resitivity of the exoplasm
  parameter SI.Resistivity exoplasm_resistivity = 300.0 * Ohm * cm;
  // Areal Capacitance of the neural membrane
  parameter SI.CapacitancePerArea membrane_capacitance_per_area = 2.0 * uF / cm2;
  // Leak conductance per unit area
  parameter Real leak_conductance_per_area = 30.365 * mS / cm2;
  // Width of the Nodes of Ranvenier
  parameter SI.Length intranodal_gap = 2.5 * um;
  // Distance between nodes as a ratio to the fiber diameter
  parameter SI.Length internodal_distance = 100 * fiber_diameter;
  // Ratio of the axon diameter to the fiber diameter
  parameter Real axon_fiber_ratio = 0.7;
  // Resting potential;
  constant SI.Voltage Vr = -70 * mV;
  // Base Sodium Permeability constant
  constant SI.Velocity basePNAB = 8e-3 * cm;
  // Base Potassium Permeability constant
  constant SI.Velocity PKB = 1.2e-3 * cm;
  // Base non-specific permeability constant
  constant SI.Velocity PPB = 0.54e-3 * cm;
  // Equilibrium leak potential
  constant SI.Voltage Vl = 0.026 * mV;
  // External Sodium concentration
  constant SI.Concentration C_Na_e = 114.5 * mmol / L;
  // Internal Sodium concentration
  constant SI.Concentration C_Na_i = 13.7 * mmol / L;
  // External Potassium concentration
  constant SI.Concentration C_K_e = 2.5 * mmol / L;
  // Internal Potassium concentration
  constant SI.Concentration C_K_i = 120.0 * mmol / L;
  parameter SI.Voltage ref_voltage = 0.0;
  parameter SI.Voltage Vth = 80 * mV;
  parameter Boolean simulate_body = false;
  parameter Integer Nnodes = 51;
  SI.Length axon_diameter = fiber_diameter * axon_fiber_ratio;
  SI.Area node_area = pi * intranodal_gap * axon_diameter;
  SI.Capacitance Cm = membrane_capacitance_per_area * node_area;
  SI.Conductance Gm = leak_conductance_per_area * node_area;
  SI.Conductance Ga = 1000.0 * pi * axon_diameter ^ 2 / (4 * axoplasm_resistivity * internodal_distance);
  // SENN has 1000 mult that makes little sense
  Real Gl = leak_conductance_per_area * pi * axon_diameter * intranodal_gap;
  Real FRT = F / R / T;
  annotation(Diagram(coordinateSystem(extent = {{-150, -90}, {150, 90}}, preserveAspectRatio = true, initialScale = 0.1, grid = {5, 5})), Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Text(visible = true, textColor = {64, 64, 64}, extent = {{-150, 60}, {150, 100}}, textString = "%name"), Rectangle(visible = true, origin = {0, -25}, lineColor = {0, 114, 195}, fillColor = {128, 202, 255}, fillPattern = FillPattern.Solid, extent = {{-100, -75}, {100, 75}}, radius = 10), Line(visible = true, points = {{-100, 0}, {100, 0}}, color = {0, 127, 255}), Line(visible = true, origin = {0, -50}, points = {{-100, 0}, {100, 0}}, color = {0, 127, 255}), Line(visible = true, origin = {0, -25}, points = {{0, 75}, {0, -75}}, color = {0, 127, 255})}));
end NeuronModelParameters;
