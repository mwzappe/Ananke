model Neuron "Model of a neuron"
  import pi = Modelica.Constants.pi;
  import SI = Modelica.SIunits;
  import R = Modelica.Constants.R;
  import F = Modelica.Constants.F;
  constant SI.Length cm = 0.01;
  constant SI.Area cm2 = SI.Conversions.from_cm2(1);
  constant SI.Capacitance uF = 1e-6;
  constant SI.Resistance Ohm = 1e-6;
  constant SI.Conductance mS = 1e-3;
  constant SI.Temperature T = 295.16;
  constant SI.Voltage Vr = -0.070;
  // Resting potential;
  constant SI.Voltage Vl = 0.0000260430075;
  // Equilibrium leak potential
  constant SI.Conductance Gl = 30.3;
  /* what units */
  constant SI.Concentration C_Na_e = 114.5;
  /* what units */
  constant SI.Concentration C_Na_i = 114.5;
  /* what units */
  constant SI.Concentration C_K_e = 114.5;
  /* what units */
  constant SI.Concentration C_K_i = 114.5;
  /* what units */
  /*   
  static constexpr senn_real_t SODO=114.5; // !(Na)o external sodium concentration
  static constexpr senn_real_t SODI=13.74; // !(Na)i internal sodium concentration
  static constexpr senn_real_t POTO=2.5; //   !(  K)o, external potassium concentration
  static constexpr senn_real_t POTI=120.; //  !(K)i, internal potassium concentration


  static constexpr senn_real_t PKB=.0012; //  !potassium permeability constant
  static constexpr senn_real_t PPB=.00054; // !nonspecific permeability constant

  static constexpr senn_real_t GL=30.3; //    !leak conductance


  static constexpr senn_real_t ELD=100.0; //
*/
  inner parameter SI.CapacitancePerArea membrane_capacitance_per_area = 2.0 * uF / cm2;
  /*F/m^2*/
  inner parameter SI.Resistivity axoplasm_resistivity = 110.0 * Ohm * cm;
  inner parameter SI.Resistivity exoplasm_resistivity = 300.0 * Ohm * cm;
  inner parameter SI.Voltage ref_voltage = 0.0;
  inner parameter SI.Voltage Vth = 0.080;
  inner parameter Real linear_membrane_conductance = 30.365 * mS / cm2;
  parameter Boolean simulate_body = false;

  model Body "Body of a neuron"
    parameter SI.Length body_width = from_cm(0.002);
    parameter SI.Length diameter = from_cm(0.002);
  end Body;

  model Hillock "Hillock of a neuron"
    parameter SI.Length diameter = from_cm(0.002);
  end Hillock;

  model Axon "Axon of a neuron"
    outer parameter SI.CapacitancePerArea membrane_capacitance_per_area;
    outer parameter SI.Resistivity axoplasm_resistivity;
    outer parameter SI.Resistivity exoplasm_resistivity;
    outer parameter SI.Voltage ref_voltage;
    outer parameter SI.Voltage Vth;
    outer parameter Real linear_membrane_conductance;
    parameter Integer Nnodes = 11;
    parameter SI.Length diameter = 0.002 / 100;
    /* cm */
    /* cm->m */
    parameter Real axon_fiber_ratio = 0.7;
    parameter SI.Length intranodal_gap = 0.00025 / 100;
    /* cm */
    /* cm->m */
    parameter SI.Length internodal_distance = 100 * diameter;
    inner SI.Area node_area = pi * intranodal_gap * diameter;
    inner SI.Capacitance Cm = membrane_capacitance_per_area * node_area;
    inner SI.Conductance Gm = linear_membrane_conductance * node_area;
    inner SI.Conductance Ga = pi * diameter ^ 2 / (4 * axoplasm_resistivity * internodal_distance);
    /*
connector Pin "Pin of an electrical component"
  SI.ElectricPotential v "Potential at the pin" annotation(unassignedMessage = "An electrical potential cannot be uniquely calculated.
The reason could be that
- a ground object is missing (Modelica.Electrical.Analog.Basic.Ground)
  to define the zero potential of the electrical circuit, or
- a connector of an electrical component is not connected.");
  flow SI.Current i "Current flowing into the pin" annotation(unassignedMessage = "An electrical current cannot be uniquely calculated.
The reason could be that
- a ground object is missing (Modelica.Electrical.Analog.Basic.Ground)
  to define the zero potential of the electrical circuit, or
- a connector of an electrical component is not connected.");
  annotation(defaultComponentName = "pin", Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}), graphics = {Rectangle(lineColor = {10, 90, 224}, fillColor = {10, 90, 224}, fillPattern = FillPattern.Solid, extent = {{-100, -100}, {100, 100}})}), Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}), graphics = {Rectangle(lineColor = {10, 90, 224}, fillColor = {10, 90, 224}, fillPattern = FillPattern.Solid, extent = {{-40, -40}, {40, 40}}), Text(origin = {0, 16.667}, textColor = {64, 64, 64}, extent = {{-160, 33.333}, {40, 73.333}}, textString = "%name")}), Documentation(revisions = "<html>
<ul>
<li><em> 1998   </em>
       by Christoph Clauss<br> initially implemented<br>
       </li>
</ul>
</html>", info = "<html>
<p>Pin is the basic electric connector. It includes the voltage which consists between the pin and the ground node. The ground node is the node of (any) ground device (Modelica.Electrical.Basic.Ground). Furthermore, the pin includes the current, which is considered to be <strong>positive</strong> if it is flowing at the pin <strong>into the device</strong>.</p>
</html>"));
end Pin;
*/

    connector NeuronSectionConn
      SI.Voltage v;
      flow SI.Current i;
    end NeuronSectionConn;

    model Axoplasm "Axoplasm between nodes"
      outer parameter SI.Conductance Ga;
      NeuronSectionConn left, right;
      SI.Voltage v;
    equation
      v = right.v - left.v;
      left.i = Ga * v;
      0 = left.i + right.i;
    end Axoplasm;

    model Node "Node of Ranvenier"
      //input SI.ElectricFieldStrength Ex;
      outer parameter SI.Conductance Gm;
      outer parameter SI.Capacitance Cm;
      input SI.Voltage Ve;
      NeuronSectionConn left, right;
      SI.Voltage V;
      SI.Current I_ionic;
      /*
      SI.Current I_Na;
      SI.Current I_K;
      SI.Current I_L;
      SI.Current I_P;
*/
    equation
      left.v = V + Ve;
      right.v = V + Ve;
      I_ionic = Gm * V;
      der(V) = (right.i - left.i) / Cm;
    end Node;

    Axoplasm a[Nnodes - 1];
    Node n[Nnodes];
  equation
    for i in 1:Nnodes loop
      n[i].Ve = 0;
    end for;
    for i in 1:Nnodes - 1 loop
      connect(n[i].right, a[i].left);
      connect(a[i].right, n[i + 1].left);
    end for;
  end Axon;

  Axon axon;
  annotation(Diagram(coordinateSystem(extent = {{-150, -90}, {150, 90}}, preserveAspectRatio = true, initialScale = 0.1, grid = {5, 5})), Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(visible = true, lineColor = {0, 114, 195}, fillColor = {255, 255, 255}, extent = {{-100, -100}, {100, 100}}, radius = 25), Text(visible = true, textColor = {64, 64, 64}, extent = {{-150, 110}, {150, 150}}, textString = "%name")}));
end Neuron;
