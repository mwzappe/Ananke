model Neuron "Model of a neuron"
  import pi = Modelica.Constants.pi;
  import SI = Modelica.SIunits;
  parameter NeuronModelParameters params;
  input SI.Voltage Ve[params.Nnodes];

  model Body "Body of a neuron"
    parameter SI.Length body_width = from_cm(0.002);
    parameter SI.Length diameter = from_cm(0.002);
  end Body;

  model Hillock "Hillock of a neuron"
    parameter SI.Length diameter = from_cm(0.002);
  end Hillock;

  model Axon "Axon of a neuron"
    parameter NeuronModelParameters params;
    input SI.Voltage Ve[params.Nnodes];

    model Node "Node of Ranvenier"
      import Modelica.Blocks.Interfaces.RealInput;
      import Modelica.Blocks.Interfaces.RealOutput;
      //input SI.ElectricFieldStrength Ex;
      parameter NeuronModelParameters params;
      input SI.Voltage Ve;
      RealInput left "left Vtot";
      RealInput right "right Vtot";
      RealOutput Vtot "Total of Trans-membrane potential and external potential";
      SI.Voltage V;
      SI.Current linear_ionic;
      Real PNAB = params.leak_conductance_per_area / params.base_leak_conductance_per_area * params.basePNAB;
      Real EFRT;

      model perm_param_form_1
        input Real a0, a1, a2;
        input String name;
        output Real value;
      equation
        assert(value > (-100) and value < 100, "Parameter " + name + " out of bounds (" + String(a0) + ", " + String(a1) + ", " + String(a2) + ")", AssertionLevel.warning);
        value = if abs(a1 / a2) < 1e-9 then a0 * a2 else a0 * a1 / (1 - exp(-a1 / a2));
      end perm_param_form_1;

      model perm_param_form_2
        input Real a0, a1, a2;
        input String name;
        output Real value;
      equation
        assert(value > (-100) and value < 100, "Parameter " + name + " out of bounds", AssertionLevel.warning);
        value = if abs(a1 / a2) < 1e-9 then a0 * a2 else a0 / (1 + exp(a1 / a2));
      end perm_param_form_2;

      model perm_param
        output Real value;
        input Real a0, a1, a2, b0, b1, b2;
        input String name;
        perm_param_form_1 alpha(a0 = a0, a1 = a1, a2 = a2, name = name + " alpha");
        perm_param_form_1 beta(a0 = b0, a1 = b1, a2 = b2, name = name + " beta");
      initial equation
        der(value) = 0;
      equation
        assert(value > (-100) and value < 100, "Parameter " + name + " out of bounds", AssertionLevel.warning);
        der(value) = 1000.0 * (alpha.value * (1 - value) - beta.value * value);
      end perm_param;

      model perm_param2
        output Real value;
        input Real a0, a1, a2, b0, b1, b2;
        input String name;
        perm_param_form_1 alpha(a0 = a0, a1 = a1, a2 = a2, name = name + " alpha");
        perm_param_form_2 beta(a0 = b0, a1 = b1, a2 = b2, name = name + " beta");
      initial equation
        der(value) = 0;
      equation
        assert(value > (-100) and value < 100, "Parameter " + name + " out of bounds", AssertionLevel.warning);
        der(value) = 1000.0 * (alpha.value * (1 - value) - beta.value * value);
      end perm_param2;

      perm_param m(name = "m", a0 = 0.36, a1 = 1000 * V - 22, a2 = 3, b0 = 0.4, b1 = 13 - 1000 * V, b2 = 20);
      perm_param2 h(name = "h", a0 = 0.1, a1 = (-10) - 1000 * V, a2 = 6, b0 = 4.5, b1 = 45 - 1000 * V, b2 = 10);
      perm_param n(name = "n", a0 = 0.02, a1 = 1000 * V - 35, a2 = 10, b0 = 0.05, b1 = 10 - 1000 * V, b2 = 10);
      perm_param p(name = "p", a0 = 0.006, a1 = 1000 * V - 40, a2 = 10, b0 = 0.09, b1 = (-25) - 1000 * V, b2 = 20);
      SI.CurrentDensity I_Na;
      SI.CurrentDensity I_K;
      SI.CurrentDensity I_L;
      SI.CurrentDensity I_P;
      Real ex;
      Real b;
      SI.Current total_ionic;
      SI.Current tim;
    initial equation
      V = 0;
    equation
      assert(der(V) > (-1e18) and der(V) < 1e18, "Trans-membrane potential derivative out-of-bounds: (" + String(der(V)) + ", " + String(tim / params.Cm) + ", " + String(params.node_area * total_ionic / params.Cm) + ")", AssertionLevel.warning);
      assert(V > (-1e5) and V < 1e5, "Trans-membrane potential out-of-bounds: " + String(V), AssertionLevel.warning);
      assert(total_ionic > (-1e5) and total_ionic < 1e5, "Total ionic current out-of-bounds: (" + String(I_Na) + ", " + String(I_K) + ", " + String(I_P) + ", " + String(I_L) + ")", AssertionLevel.warning);
      assert(tim > (-1e5) and tim < 1e5, "Total axoplasmic current out-of-bounds: (" + String(left) + " | " + String(V) + ", " + String(Ve) + " | " + String(right) + ")", AssertionLevel.error);
      Vtot = V + Ve;
      linear_ionic = params.Gm * V;
      // Linear model
      // FRT->mV-1
      EFRT = (V + params.Vr) * params.FRT;
      // Unitless
      ex = exp(EFRT);
      // SENN -> C mmol-1
      b = EFRT * params.F / (1 - ex);
      // SENN -> 1000.0 * (cm s-1) (C mmol-1) (mmol L-1) *NOTE* L wasn't accounted for so this is 1000x
      // So add another 1000.0 to each of these to account for the fix to the Liter unit
      // 1000.0 * 1000.0 * (cm s-1) (C mmol-1) (mmol L-1) -> 1e6 (cm s-1) (C 1e-3 cm-3) -> 1e6 (C s-1) (1e-3 cm-2) -> uA / (1e-3 cm-2)
      // However, this *still* doesn't quite work right, it should be uA cm-2
      // But the permeability constant is given in cm s-1, when the time unit is ms, so....
      // 1e6 (C 1000 ms) -> 1e9/1e3->1e6, ergo uA/cm-2
      I_Na = PNAB * h.value * m.value ^ 2 * b * (params.C_Na_e - params.C_Na_i * ex);
      I_K = params.PKB * n.value ^ 2 * b * (params.C_K_e - params.C_K_i * ex);
      I_P = params.PPB * p.value ^ 2 * b * (params.C_Na_e - params.C_Na_i * ex);
      // SENN -> mS cm-2 mV -> uA cm-2
      I_L = params.leak_conductance_per_area * (V - params.Vl);
      // SENN -> mA cm-2
      total_ionic = params.node_area * (I_Na + I_K + I_P + I_L);
      // SENN -> mS mV -> uA
      tim = params.Ga * (right - 2 * Vtot + left);
      // SENN -> (mA - cm2 * (mA cm-2)) / uF -> 1e6 V/s
      der(V) = (tim - total_ionic) / params.Cm;
    end Node;

    SI.Voltage center_diff;
    SI.Current center_current;
    Node n[params.Nnodes](params = params);
  equation
    center_diff = n[5].Vtot - n[4].Vtot - (n[4].Vtot - n[3].Vtot);
    center_current = params.Ga * center_diff;
    for i in 1:params.Nnodes loop
      n[i].Ve = Ve[i];
    end for;
    for i in 1:params.Nnodes - 1 loop
      connect(n[i].right, n[i + 1].Vtot);
      connect(n[i].Vtot, n[i + 1].left);
    end for;
    n[1].left = n[1].Vtot;
    n[params.Nnodes].right = n[params.Nnodes].Vtot;
  end Axon;

  Axon axon(params = params, Ve = Ve);
  annotation(experiment(StopTime = 0.001), Diagram(coordinateSystem(extent = {{-150, -90}, {150, 90}}, preserveAspectRatio = true, initialScale = 0.1, grid = {5, 5})), Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(visible = true, lineColor = {0, 114, 195}, fillColor = {255, 255, 255}, extent = {{-100, -100}, {100, 100}}, radius = 25), Text(visible = true, textColor = {64, 64, 64}, extent = {{-150, 110}, {150, 150}}, textString = "%name")}));
end Neuron;
