model TwoElectrodeTest
  import pi = Modelica.Constants.pi;
  constant NeuronModelParameters params(Nnodes = 51);
  TwoElectrodeStimulus stim(params = params);
  Neuron neuron(params = params);
equation
  stim.I = 16 * sin(2 * pi * 5 * time);
  neuron.Ve = stim.V;
  annotation(experiment(StopTime = 6e-6, Interval = 200e-9, __Wolfram_Algorithm = "rk4"), Diagram(coordinateSystem(extent = {{-150, -90}, {150, 90}}, preserveAspectRatio = true, initialScale = 0.1, grid = {5, 5})), Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(visible = true, lineColor = {0, 114, 195}, fillColor = {255, 255, 255}, extent = {{-100, -100}, {100, 100}}, radius = 25), Text(visible = true, textColor = {64, 64, 64}, extent = {{-150, 110}, {150, 150}}, textString = "%name")}));
end TwoElectrodeTest;
