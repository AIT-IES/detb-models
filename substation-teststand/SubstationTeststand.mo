within ;
package SubstationTeststand

  model TestStandDummy
    extends Components.BaseClasses.TeststandBase(
      redeclare Components.SecondarySide.DummySourceAndSink secondarySide(
          m_flow_nominal=4),
      redeclare Components.Substation.Substation substation(m1_flow_nominal=4,
          m2_flow_nominal=4),
      redeclare Components.PrimarySide.DummySourceAndSink primarySide(
          m_flow_nominal=4));

    annotation (Icon(graphics={Rectangle(
            extent={{-100,100},{100,-100}},
            lineColor={0,0,0},
            fillColor={128,255,0},
            fillPattern=FillPattern.Solid,
            radius=10)}));
  end TestStandDummy;

  model TestStandSimplifiedDirectPressureCtrl
    extends Components.BaseClasses.TeststandBase(
      redeclare Components.SecondarySide.SimplifiedDirectPressureCtrl
        secondarySide(m_flow_nominal=4),
      redeclare Components.Substation.Substation substation(m1_flow_nominal=4,
          m2_flow_nominal=4),
      redeclare Components.PrimarySide.SimplifiedDirectPressureCtrl primarySide(
          m_flow_nominal=4));

    annotation (Icon(graphics={Rectangle(
            extent={{-100,100},{100,-100}},
            lineColor={0,0,0},
            fillColor={93,180,0},
            fillPattern=FillPattern.Solid,
            radius=10)}));
  end TestStandSimplifiedDirectPressureCtrl;

  model TestStandSimplifiedPIDPressureCtrl
    extends Components.BaseClasses.TeststandBase(
      redeclare Components.SecondarySide.SimplifiedPIDPressureCtrl
        secondarySide(m_flow_nominal=4),
      redeclare Components.Substation.Substation substation(m1_flow_nominal=4,
          m2_flow_nominal=4),
      redeclare Components.PrimarySide.SimplifiedPIDPressureCtrl primarySide(
          m_flow_nominal=4));

    annotation (Icon(graphics={Rectangle(
            extent={{-100,100},{100,-100}},
            lineColor={0,0,0},
            fillColor={51,100,0},
            fillPattern=FillPattern.Solid,
            radius=10)}));
  end TestStandSimplifiedPIDPressureCtrl;

  package Components

    package Substation
      model Substation
        extends BaseClasses.SubstationBase;

        IBPSA.Fluid.HeatExchangers.ConstantEffectiveness hex(
          redeclare package Medium1 = IBPSA.Media.Water,
          redeclare package Medium2 = IBPSA.Media.Water,
          m1_flow_nominal=m1_flow_nominal,
          m2_flow_nominal=m2_flow_nominal,
          dp1_nominal(displayUnit="bar") = dp1_nominal,
          dp2_nominal(displayUnit="bar") = dp2_nominal)
          annotation (Placement(transformation(extent={{10,-10},{30,10}})));
        IBPSA.Fluid.Actuators.Valves.TwoWayLinear val(
          redeclare package Medium = IBPSA.Media.Water,
          m_flow_nominal=m1_flow_nominal,
          dpValve_nominal(displayUnit="bar") = dpValve_nominal,
          riseTime=10,
          y_start=0,
          dpFixed_nominal(displayUnit="bar") = 0)
          annotation (Placement(transformation(extent={{-60,70},{-40,50}})));
        Modelica.Blocks.Continuous.LimPID PID(
          controllerType=Modelica.Blocks.Types.SimpleController.PID,
          k=0.003,
          yMax=1,
          yMin=0) annotation (Placement(transformation(extent={{-10,-10},{10,10}},
              rotation=90,
              origin={-50,30})));
        IBPSA.Fluid.Sensors.Temperature T_supply_secondary_sense_K(redeclare
            package
            Medium = IBPSA.Media.Water) annotation (Placement(transformation(
              extent={{10,-10},{-10,10}},
              rotation=0,
              origin={-20,30})));


      equation
        connect(hex.port_b1, port_b1) annotation (Line(points={{30,6},{40,6},{40,60},{
                100,60}},      color={0,127,255}));
        connect(port_a2, hex.port_a2) annotation (Line(points={{100,-60},{40,-60},{40,
                -6},{30,-6}},     color={0,127,255}));
        connect(port_b2, hex.port_b2) annotation (Line(points={{-100,-60},{0,-60},{0,-6},
                {10,-6}},                color={0,127,255}));
        connect(T_supply_secondary_sense_K.port, hex.port_b2)
          annotation (Line(points={{-20,20},{-20,-6},{10,-6}}, color={0,127,255}));
        connect(port_a1, val.port_a)
          annotation (Line(points={{-100,60},{-60,60}}, color={0,127,255}));
        connect(val.port_b, hex.port_a1) annotation (Line(points={{-40,60},{0,60},{0,6},
                {10,6}}, color={0,127,255}));
        connect(val.y, PID.y) annotation (Line(points={{-50,48},{-50,41}},
                               color={0,0,127}));
        connect(PID.u_m, T_supply_secondary_sense_K.T)
          annotation (Line(points={{-38,30},{-27,30}}, color={0,0,127}));
        connect(T_supply_secondary_set, PID.u_s) annotation (Line(points={{-120,
                0},{-50,0},{-50,18}}, color={0,0,127}));
        annotation (Icon(graphics={Rectangle(
                extent={{-100,100},{100,-100}},
                lineColor={0,0,0},
                fillColor={28,108,200},
                fillPattern=FillPattern.Solid,
                radius=10)}));
      end Substation;
    end Substation;

    package PrimarySide
      model DummySourceAndSink
        extends BaseClasses.PrimarySideBase;
        IBPSA.Fluid.Sources.Boundary_pT sink_primary(
          redeclare package Medium = IBPSA.Media.Water,
          nPorts=1) annotation (Placement(transformation(extent={{-40,-10},{-60,
                  10}})));
        IBPSA.Fluid.Sources.Boundary_pT source_primary(
          redeclare package Medium = IBPSA.Media.Water,
          use_p_in=true,
          use_T_in=true,
          nPorts=1)
          annotation (Placement(transformation(extent={{40,-10},{60,10}})));
        Modelica.Blocks.Sources.RealExpression p_primary_supply_Pa(y=
              delta_p_primary_set + IBPSA.Media.Water.p_default)
          annotation (Placement(transformation(extent={{-20,-70},{0,-50}})));
      equation
        connect(port_a, sink_primary.ports[1])
          annotation (Line(points={{-100,0},{-60,0}}, color={0,127,255}));
        connect(source_primary.ports[1], port_b)
          annotation (Line(points={{60,0},{100,0}}, color={0,127,255}));
        connect(source_primary.p_in, p_primary_supply_Pa.y) annotation (Line(
              points={{38,8},{18,8},{18,-60},{1,-60}},   color={0,0,127}));
        connect(T_supply_primary_set, source_primary.T_in) annotation (Line(
              points={{40,-120},{40,-60},{22,-60},{22,4},{38,4}},
                                                         color={0,0,127}));
        annotation (Icon(graphics={Rectangle(
                extent={{-100,100},{100,-100}},
                lineColor={0,0,0},
                fillColor={128,255,0},
                fillPattern=FillPattern.Solid,
                radius=10)}));
      end DummySourceAndSink;

      model SimplifiedDirectPressureCtrl
        extends BaseClasses.PrimarySideBase;
        IBPSA.Fluid.HeatExchangers.Heater_T hea(
          redeclare package Medium = IBPSA.Media.Water,
          m_flow_nominal=m_flow_nominal,
          dp_nominal=2000,
          T_start=353.15,
          energyDynamics=Modelica.Fluid.Types.Dynamics.SteadyState)
          annotation (Placement(transformation(extent={{10,10},{30,-10}})));
        IBPSA.Fluid.Movers.FlowControlled_dp pump(
          redeclare package Medium = IBPSA.Media.Water,
          T_start=353.15,
          m_flow_nominal=m_flow_nominal,
          redeclare IBPSA.Fluid.Movers.Data.Generic per,
          nominalValuesDefineDefaultPressureCurve=true,
          dp_nominal(displayUnit="bar") = 400000)
          annotation (Placement(transformation(extent={{-30,10},{-10,-10}})));
        IBPSA.Fluid.Storage.ExpansionVessel exp(redeclare package Medium =
              IBPSA.Media.Water, V_start=1)
          annotation (Placement(transformation(extent={{-50,20},{-30,40}})));
        IBPSA.Fluid.FixedResistances.PlugFlowPipe pipe_return(
          redeclare package Medium = IBPSA.Media.Water,
          dh=0.03,
          v_nominal=4,
          length=5,
          m_flow_nominal=m_flow_nominal,
          dIns=0.01,
          kIns=0.02,
          T_start_in=313.15,
          nPorts=1) annotation (Placement(transformation(
              extent={{10,10},{-10,-10}},
              rotation=180,
              origin={-70,0})));
        Modelica.Thermal.HeatTransfer.Sources.FixedTemperature
          ambientTemperature(T=293.15) annotation (Placement(transformation(
              extent={{10,-10},{-10,10}},
              rotation=90,
              origin={0,70})));
        IBPSA.Fluid.FixedResistances.PlugFlowPipe pipe_supply(
          redeclare package Medium = IBPSA.Media.Water,
          dh=0.03,
          v_nominal=4,
          length=5,
          m_flow_nominal=m_flow_nominal,
          dIns=0.01,
          kIns=0.02,
          T_start_in=353.15,
          nPorts=1) annotation (Placement(transformation(
              extent={{10,10},{-10,-10}},
              rotation=180,
              origin={70,0})));
      equation
        connect(pump.port_b, hea.port_a)
          annotation (Line(points={{-10,0},{10,0}}, color={0,127,255}));
        connect(pump.port_a, exp.port_a) annotation (Line(points={{-30,0},{-40,
                0},{-40,20}}, color={0,127,255}));
        connect(hea.port_b, pipe_supply.port_a) annotation (Line(points={{30,0},
                {46,0},{46,1.77636e-15},{60,1.77636e-15}}, color={0,127,255}));
        connect(pipe_supply.ports_b[1], port_b) annotation (Line(points={{80,
                -4.44089e-16},{89,-4.44089e-16},{89,0},{100,0}}, color={0,127,
                255}));
        connect(port_a, pipe_return.port_a) annotation (Line(points={{-100,0},{
                -90,0},{-90,1.77636e-15},{-80,1.77636e-15}}, color={0,127,255}));
        connect(pipe_return.ports_b[1], pump.port_a) annotation (Line(points={{
                -60,-4.44089e-16},{-46,-4.44089e-16},{-46,0},{-30,0}}, color={0,
                127,255}));
        connect(pipe_return.heatPort, ambientTemperature.port) annotation (Line(
              points={{-70,10},{-70,50},{-4.44089e-16,50},{-4.44089e-16,60}},
              color={191,0,0}));
        connect(pipe_supply.heatPort, ambientTemperature.port) annotation (Line(
              points={{70,10},{70,50},{-4.44089e-16,50},{-4.44089e-16,60}},
              color={191,0,0}));
        connect(delta_p_primary_set, pump.dp_in) annotation (Line(points={{-40,
                -120},{-40,-20},{-20,-20},{-20,-12}},
                                           color={0,0,127}));
        connect(T_supply_primary_set, hea.TSet) annotation (Line(points={{40,-120},
                {40,-20},{0,-20},{0,-8},{8,-8}},
                                             color={0,0,127}));
        annotation (Icon(graphics={Rectangle(
                extent={{-100,100},{100,-100}},
                lineColor={0,0,0},
                fillColor={93,180,0},
                fillPattern=FillPattern.Solid,
                radius=10)}));
      end SimplifiedDirectPressureCtrl;

      model SimplifiedPIDPressureCtrl
        extends BaseClasses.PrimarySideBase;

        parameter Modelica.SIunits.PressureDifference pump_dp_nominal = 1000000
          "Nominal pressure difference of pump at full speed";

        parameter Real pump_m_flow_max=5
          "Maximum mass flow generated by pump at full speed";

        parameter Real pump_speed_y_set = 0.5
          "Relative speed of pump (a value of 1 will generate maximum mass flow)";

        IBPSA.Fluid.HeatExchangers.Heater_T hea(
          redeclare package Medium = IBPSA.Media.Water,
          m_flow_nominal=m_flow_nominal,
          dp_nominal=2000,
          T_start=353.15)
          annotation (Placement(transformation(extent={{50,60},{70,40}})));
        IBPSA.Fluid.Movers.FlowControlled_m_flow
                                             pump(
          redeclare package Medium = IBPSA.Media.Water,
          T_start=313.15,
          m_flow_nominal=pump_m_flow_max,
          m_flow_small=0.01,
          per(pressure(V_flow={0,pump_m_flow_max}/1.2,
                               dp={pump_dp_nominal,0})),
          riseTime=10,
          dp_nominal=pump_dp_nominal)
          annotation (Placement(transformation(extent={{-8,42},{8,58}})));
        IBPSA.Fluid.Storage.ExpansionVessel exp(redeclare package Medium =
              IBPSA.Media.Water, V_start=1)
          annotation (Placement(transformation(extent={{-28,62},{-12,78}})));
        IBPSA.Fluid.FixedResistances.Junction jun_p_low(
          redeclare package Medium = IBPSA.Media.Water,
          T_start=313.15,
          m_flow_nominal={m_flow_nominal,-m_flow_nominal,-m_flow_nominal},
          dp_nominal(displayUnit="Pa") = {1000,-1000,-1000})
          annotation (Placement(transformation(extent={{-56,44},{-44,56}})));
        IBPSA.Fluid.FixedResistances.Junction jun_p_high(
          redeclare package Medium = IBPSA.Media.Water,
          T_start=313.15,
          m_flow_nominal={m_flow_nominal,-m_flow_nominal,m_flow_nominal},
          dp_nominal={1000,-1000,1000})
          annotation (Placement(transformation(extent={{24,44},{36,56}})));
        IBPSA.Fluid.Actuators.Valves.TwoWayLinear valve(
          redeclare package Medium = IBPSA.Media.Water,
          m_flow_nominal=m_flow_nominal,
          dpValve_nominal(displayUnit="bar") = 400000,
          riseTime=10,
          y_start=0,
          dpFixed_nominal(displayUnit="bar") = 0)
          annotation (Placement(transformation(extent={{-18,8},{-2,-8}})));
        IBPSA.Fluid.Sensors.RelativePressure senRelPre(redeclare package Medium =
              IBPSA.Media.Water)
          annotation (Placement(transformation(extent={{-2,28},{-18,12}})));
        Modelica.Blocks.Sources.RealExpression pump_speed_y(y=pump_speed_y_set*
              pump_m_flow_max)
          annotation (Placement(transformation(
              extent={{10,-10},{-10,10}},
              rotation=0,
              origin={50,80})));
        IBPSA.Fluid.FixedResistances.PlugFlowPipe pipe_return(
          redeclare package Medium = IBPSA.Media.Water,
          dh=0.03,
          v_nominal=4,
          length=5,
          m_flow_nominal=m_flow_nominal,
          dIns=0.01,
          kIns=0.02,
          T_start_in=313.15,
          nPorts=1) annotation (Placement(transformation(
              extent={{-8,8},{8,-8}},
              rotation=90,
              origin={-80,20})));
        IBPSA.Fluid.FixedResistances.PlugFlowPipe pipe_supply(
          redeclare package Medium = IBPSA.Media.Water,
          dh=0.03,
          v_nominal=4,
          length=5,
          m_flow_nominal=m_flow_nominal,
          dIns=0.01,
          kIns=0.02,
          T_start_in=353.15,
          nPorts=1) annotation (Placement(transformation(
              extent={{8,-8},{-8,8}},
              rotation=90,
              origin={80,20})));
        Modelica.Thermal.HeatTransfer.Sources.FixedTemperature
          ambientTemperature(T=293.15) annotation (Placement(transformation(
              extent={{-8,-8},{8,8}},
              rotation=90,
              origin={10,-50})));

        Modelica.Blocks.Sources.RealExpression delta_p_primary_sense_Pa(y=senRelPre.p_rel)
          annotation (Placement(transformation(extent={{-90,-80},{-70,-60}})));
        IBPSA.Controls.Continuous.LimPID PID_delta_p(
          controllerType=Modelica.Blocks.Types.SimpleController.PID,
          k=8e-6,
          Ti=1,
          yMax=1,
          yMin=0,
          reverseAction=true) annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=90,
              origin={-40,-70})));

      initial equation
        assert(pump_speed_y_set >= 0 and pump_speed_y_set <= 1,
         "Relative speed of pump must be between 0 and 1");

      equation
        connect(jun_p_low.port_2, pump.port_a)
          annotation (Line(points={{-44,50},{-8,50}}, color={0,127,255}));
        connect(pump.port_b, jun_p_high.port_1)
          annotation (Line(points={{8,50},{24,50}}, color={0,127,255}));
        connect(jun_p_high.port_2, hea.port_a)
          annotation (Line(points={{36,50},{50,50}}, color={0,127,255}));
        connect(jun_p_low.port_3, valve.port_a) annotation (Line(points={{-50,44},{-50,
                0},{-18,0}},          color={0,127,255}));
        connect(valve.port_b, jun_p_high.port_3)
          annotation (Line(points={{-2,0},{30,0},{30,44}}, color={0,127,255}));
        connect(port_a, pipe_return.port_a) annotation (Line(points={{-100,0},{
                -80,0},{-80,12}}, color={0,127,255}));
        connect(pipe_return.ports_b[1], jun_p_low.port_1) annotation (Line(
              points={{-80,28},{-80,50},{-56,50}}, color={0,127,255}));
        connect(port_b, pipe_supply.ports_b[1]) annotation (Line(points={{100,0},
                {80,0},{80,12}}, color={0,127,255}));
        connect(pipe_supply.port_a, hea.port_b) annotation (Line(points={{80,28},{80,50},
                {70,50}},         color={0,127,255}));
        connect(pipe_return.heatPort, ambientTemperature.port) annotation (Line(
              points={{-72,20},{-60,20},{-60,-20},{10,-20},{10,-42}}, color={
                191,0,0}));
        connect(pipe_supply.heatPort, ambientTemperature.port) annotation (Line(
              points={{72,20},{60,20},{60,-20},{10,-20},{10,-42}}, color={191,0,
                0}));
        connect(exp.port_a, pump.port_a) annotation (Line(points={{-20,62},{-20,50},{-8,
                50}},         color={0,127,255}));
        connect(senRelPre.port_b, jun_p_low.port_2) annotation (Line(points={{-18,20},
                {-40,20},{-40,50},{-44,50}},         color={0,127,255}));
        connect(senRelPre.port_a, jun_p_high.port_1) annotation (Line(points={{-2,20},
                {20,20},{20,50},{24,50}},       color={0,127,255}));
        connect(pump_speed_y.y, pump.m_flow_in)
          annotation (Line(points={{39,80},{0,80},{0,59.6}}, color={0,0,127}));
        connect(delta_p_primary_sense_Pa.y, PID_delta_p.u_m)
          annotation (Line(points={{-69,-70},{-52,-70}}, color={0,0,127}));
        connect(delta_p_primary_set, PID_delta_p.u_s)
          annotation (Line(points={{-40,-120},{-40,-82}}, color={0,0,127}));
        connect(T_supply_primary_set, hea.TSet) annotation (Line(points={{40,-120},
                {40,-120},{40,42},{48,42}},     color={0,0,127}));
        connect(PID_delta_p.y, valve.y) annotation (Line(points={{-40,-59},{-40,
                -40},{-10,-40},{-10,-9.6}}, color={0,0,127}));
        annotation (Icon(graphics={Rectangle(
                extent={{-100,100},{100,-100}},
                lineColor={0,0,0},
                fillColor={51,100,0},
                fillPattern=FillPattern.Solid,
                radius=10)}));
      end SimplifiedPIDPressureCtrl;
    end PrimarySide;

    package SecondarySide
      model DummySourceAndSink
        extends BaseClasses.SecondarySideBase;

        IBPSA.Fluid.Sources.Boundary_pT sink_secondary(redeclare package Medium =
              IBPSA.Media.Water, nPorts=1)
                    annotation (Placement(transformation(extent={{-50,-10},{-70,
                  10}})));
        IBPSA.Fluid.Sources.Boundary_pT source_secondary(
          redeclare package Medium = IBPSA.Media.Water,
          use_p_in=true,
          use_T_in=true,
          nPorts=1)
          annotation (Placement(transformation(extent={{10,-10},{30,10}})));
        IBPSA.Fluid.Sensors.MassFlowRate senMasFlo(redeclare package Medium =
              IBPSA.Media.Water)
          annotation (Placement(transformation(extent={{40,-10},{60,10}})));
        IBPSA.Fluid.Actuators.Valves.TwoWayLinear val(
          redeclare package Medium = IBPSA.Media.Water,
          m_flow_nominal=m_flow_nominal,
          dpValve_nominal(displayUnit="bar") = 100000,
          riseTime=10,
          y_start=0,
          dpFixed_nominal(displayUnit="bar") = 50000)
          annotation (Placement(transformation(extent={{70,-10},{90,10}})));
        Modelica.Blocks.Continuous.LimPID PID(
          controllerType=Modelica.Blocks.Types.SimpleController.PID,
          k=0.03,
          yMax=1,
          yMin=0) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=-90,
              origin={80,48})));
        Modelica.Blocks.Sources.RealExpression p_secondary_supply_Pa(y=
              delta_p_secondary_set + IBPSA.Media.Water.p_default)
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=-90,
              origin={0,50})));
      equation
        connect(senMasFlo.m_flow, PID.u_m)
          annotation (Line(points={{50,11},{50,48},{68,48}}, color={0,0,127}));
        connect(PID.y, val.y)
          annotation (Line(points={{80,37},{80,12}}, color={0,0,127}));
        connect(source_secondary.ports[1], senMasFlo.port_a)
          annotation (Line(points={{30,0},{40,0}}, color={0,127,255}));
        connect(val.port_a, senMasFlo.port_b)
          annotation (Line(points={{70,0},{60,0}}, color={0,127,255}));
        connect(port_a, sink_secondary.ports[1])
          annotation (Line(points={{-100,0},{-70,0}}, color={0,127,255}));
        connect(val.port_b, port_b)
          annotation (Line(points={{90,0},{100,0}}, color={0,127,255}));
        connect(m_flow_return_secondary_set, PID.u_s)
          annotation (Line(points={{80,120},{80,60}}, color={0,0,127}));
        connect(p_secondary_supply_Pa.y, source_secondary.p_in)
          annotation (Line(points={{0,39},{0,8},{8,8}}, color={0,0,127}));
        connect(source_secondary.T_in, T_return_secondary_set) annotation (Line(
              points={{8,4},{-4,4},{-4,20},{-80,20},{-80,120}}, color={0,0,127}));
        annotation (Icon(graphics={Rectangle(
                extent={{-100,100},{100,-100}},
                lineColor={0,0,0},
                fillColor={128,255,0},
                fillPattern=FillPattern.Solid,
                radius=10)}));
      end DummySourceAndSink;

      model SimplifiedDirectPressureCtrl
        extends BaseClasses.SecondarySideBase;

        IBPSA.Fluid.Sensors.MassFlowRate senMasFlo(redeclare package Medium =
              IBPSA.Media.Water)
          annotation (Placement(transformation(extent={{22,30},{42,50}})));
        IBPSA.Fluid.Actuators.Valves.TwoWayLinear val(
          redeclare package Medium = IBPSA.Media.Water,
          m_flow_nominal=m_flow_nominal,
          dpValve_nominal(displayUnit="bar") = 100000,
          riseTime=10,
          y_start=0,
          dpFixed_nominal(displayUnit="bar") = 50000)
          annotation (Placement(transformation(extent={{50,30},{70,50}})));
        Modelica.Blocks.Continuous.LimPID PID(
          controllerType=Modelica.Blocks.Types.SimpleController.PID,
          k=0.8,
          yMax=1,
          yMin=0) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=-90,
              origin={80,80})));
        IBPSA.Fluid.Movers.FlowControlled_dp pump(
          redeclare package Medium = IBPSA.Media.Water,
          m_flow_nominal=m_flow_nominal,
          redeclare IBPSA.Fluid.Movers.Data.Generic per,
          nominalValuesDefineDefaultPressureCurve=true,
          riseTime=10,
          dp_nominal(displayUnit="bar") = 200000)
          annotation (Placement(transformation(extent={{-10,30},{10,50}})));
        IBPSA.Fluid.Storage.ExpansionVessel exp(redeclare package Medium =
              IBPSA.Media.Water, V_start=1)
          annotation (Placement(transformation(extent={{-40,70},{-20,90}})));
        IBPSA.Fluid.FixedResistances.PlugFlowPipe pipe_return(
          redeclare package Medium = IBPSA.Media.Water,
          dh=0.03,
          v_nominal=4,
          length=5,
          m_flow_nominal=m_flow_nominal,
          dIns=0.01,
          kIns=0.02,
          nPorts=1) annotation (Placement(transformation(
              extent={{10,-10},{-10,10}},
              rotation=270,
              origin={-80,20})));
        Modelica.Thermal.HeatTransfer.Sources.FixedTemperature
          ambientTemperature(T=293.15) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=90,
              origin={0,-10})));
        IBPSA.Fluid.FixedResistances.PlugFlowPipe pipe_supply(
          redeclare package Medium = IBPSA.Media.Water,
          dh=0.03,
          v_nominal=4,
          length=5,
          m_flow_nominal=m_flow_nominal,
          dIns=0.01,
          kIns=0.02,
          nPorts=1) annotation (Placement(transformation(
              extent={{10,-10},{-10,10}},
              rotation=90,
              origin={80,20})));
        IBPSA.Fluid.HeatExchangers.SensibleCooler_T coo(
          redeclare package Medium = IBPSA.Media.Water,
          m_flow_nominal=m_flow_nominal,
          dp_nominal=2000)
          annotation (Placement(transformation(extent={{-60,30},{-40,50}})));
      equation
        connect(senMasFlo.m_flow, PID.u_m)
          annotation (Line(points={{32,51},{32,80},{68,80}}, color={0,0,127}));
        connect(PID.y, val.y)
          annotation (Line(points={{80,69},{80,60},{60,60},{60,52}},
                                                     color={0,0,127}));
        connect(val.port_a, senMasFlo.port_b)
          annotation (Line(points={{50,40},{42,40}},
                                                   color={0,127,255}));
        connect(m_flow_return_secondary_set, PID.u_s)
          annotation (Line(points={{80,120},{80,92}}, color={0,0,127}));
        connect(pump.port_a,exp. port_a) annotation (Line(points={{-10,40},{-30,
                40},{-30,70}},color={0,127,255}));
        connect(port_a, pipe_return.port_a) annotation (Line(points={{-100,0},{
                -80,0},{-80,10}}, color={0,127,255}));
        connect(val.port_b, pipe_supply.port_a) annotation (Line(points={{70,40},
                {80,40},{80,30}}, color={0,127,255}));
        connect(port_b, pipe_supply.ports_b[1]) annotation (Line(points={{100,0},
                {80,0},{80,10}}, color={0,127,255}));
        connect(senMasFlo.port_a, pump.port_b)
          annotation (Line(points={{22,40},{10,40}}, color={0,127,255}));
        connect(pipe_supply.heatPort, ambientTemperature.port) annotation (Line(
              points={{70,20},{6.66134e-16,20},{6.66134e-16,0}}, color={191,0,0}));
        connect(pipe_return.heatPort, ambientTemperature.port) annotation (Line(
              points={{-70,20},{6.66134e-16,20},{6.66134e-16,0}}, color={191,0,
                0}));
        connect(coo.port_b, pump.port_a)
          annotation (Line(points={{-40,40},{-10,40}}, color={0,127,255}));
        connect(coo.port_a, pipe_return.ports_b[1]) annotation (Line(points={{
                -60,40},{-80,40},{-80,30}}, color={0,127,255}));
        connect(T_return_secondary_set, coo.TSet) annotation (Line(points={{-80,
                120},{-80,48},{-62,48}}, color={0,0,127}));
        connect(delta_p_secondary_set, pump.dp_in)
          annotation (Line(points={{0,120},{0,52}}, color={0,0,127}));
        annotation (Icon(graphics={Rectangle(
                extent={{-100,100},{100,-100}},
                lineColor={0,0,0},
                fillColor={93,180,0},
                fillPattern=FillPattern.Solid,
                radius=10)}));
      end SimplifiedDirectPressureCtrl;

      model SimplifiedPIDPressureCtrl
        extends BaseClasses.SecondarySideBase;

        IBPSA.Fluid.Sensors.MassFlowRate senMasFlo(redeclare package Medium =
              IBPSA.Media.Water)
          annotation (Placement(transformation(extent={{44,-76},{56,-64}})));
        IBPSA.Fluid.Actuators.Valves.TwoWayLinear val_m_flow(
          redeclare package Medium = IBPSA.Media.Water,
          m_flow_nominal=m_flow_nominal,
          dpValve_nominal(displayUnit="bar") = 100000,
          riseTime=10,
          y_start=0,
          dpFixed_nominal(displayUnit="bar") = 50000)
          annotation (Placement(transformation(extent={{64,-76},{76,-64}})));
        IBPSA.Controls.Continuous.LimPID PID_m_flow(
          controllerType=Modelica.Blocks.Types.SimpleController.PID,
          k=0.05,
          yMax=1,
          yMin=0) annotation (Placement(transformation(
              extent={{-6,-6},{6,6}},
              rotation=-90,
              origin={80,50})));
        IBPSA.Fluid.Movers.FlowControlled_m_flow pump(
          redeclare package Medium = IBPSA.Media.Water,
          T_start=333.15,
          m_flow_nominal=pump_m_flow_max,
          per(pressure(V_flow={0,pump_m_flow_max}/1.2,
                               dp={pump_dp_nominal,0})),
          addPowerToMedium=true,
          riseTime=10,
          dp_nominal=pump_dp_nominal)
          annotation (Placement(transformation(extent={{-26,-64},{-14,-76}})));
        IBPSA.Fluid.Storage.ExpansionVessel exp(redeclare package Medium =
              IBPSA.Media.Water, V_start=1)
          annotation (Placement(transformation(extent={{-46,-56},{-34,-44}})));
        IBPSA.Fluid.FixedResistances.PlugFlowPipe pipe_return(
          redeclare package Medium = IBPSA.Media.Water,
          dh=0.03,
          v_nominal=4,
          length=5,
          m_flow_nominal=m_flow_nominal,
          dIns=0.01,
          kIns=0.02,
          T_start_in=333.15,
          nPorts=1) annotation (Placement(transformation(
              extent={{-6,-6},{6,6}},
              rotation=270,
              origin={-80,-20})));
        Modelica.Thermal.HeatTransfer.Sources.FixedTemperature ambientTemperature(T=298.15)
          annotation (Placement(transformation(
              extent={{6,-6},{-6,6}},
              rotation=90,
              origin={-50,10})));
        IBPSA.Fluid.FixedResistances.PlugFlowPipe pipe_supply(
          redeclare package Medium = IBPSA.Media.Water,
          dh=0.03,
          v_nominal=4,
          length=5,
          m_flow_nominal=m_flow_nominal,
          dIns=0.01,
          kIns=0.02,
          T_start_in=313.15,
          nPorts=1) annotation (Placement(transformation(
              extent={{-6,-6},{6,6}},
              rotation=90,
              origin={80,-20})));
        IBPSA.Fluid.HeatExchangers.SensibleCooler_T coo(
          redeclare package Medium = IBPSA.Media.Water,
          m_flow_nominal=m_flow_nominal,
          dp_nominal=2000,
          T_start=313.15)
          annotation (Placement(transformation(extent={{24,-76},{36,-64}})));
        IBPSA.Fluid.Sensors.RelativePressure senRelPre(redeclare package Medium =
              IBPSA.Media.Water)
          annotation (Placement(transformation(extent={{-24,-24},{-36,-36}})));
        IBPSA.Fluid.FixedResistances.Junction jun_p_low(
          redeclare package Medium = IBPSA.Media.Water,
          T_start=333.15,
          m_flow_nominal={m_flow_nominal,-m_flow_nominal,-m_flow_nominal},
          dp_nominal(displayUnit="Pa") = {1000,-1000,-1000})
          annotation (Placement(transformation(extent={{-74,-66},{-66,-74}})));
        IBPSA.Fluid.FixedResistances.Junction jun_p_high(
          redeclare package Medium = IBPSA.Media.Water,
          T_start=333.15,
          m_flow_nominal={m_flow_nominal,-m_flow_nominal,m_flow_nominal},
          dp_nominal={1000,-1000,1000})
          annotation (Placement(transformation(extent={{6,-66},{14,-74}})));
        IBPSA.Controls.Continuous.LimPID PID_delta_p(
          controllerType=Modelica.Blocks.Types.SimpleController.PID,
          k=5e-6,
          Ti=1,
          yMax=1,
          yMin=0,
          reverseAction=true) annotation (Placement(transformation(
              extent={{-6,-6},{6,6}},
              rotation=-90,
              origin={0,50})));
        IBPSA.Fluid.Actuators.Valves.TwoWayLinear val_delta_p(
          redeclare package Medium = IBPSA.Media.Water,
          m_flow_nominal=m_flow_nominal,
          dpValve_nominal(displayUnit="bar") = 400000,
          riseTime=10,
          y_start=0,
          dpFixed_nominal(displayUnit="bar") = 0)
          annotation (Placement(transformation(extent={{-16,-16},{-4,-4}})));
        Modelica.Blocks.Sources.RealExpression pump_speed_y(y=pump_speed_y_set*
              pump_m_flow_max)
          annotation (Placement(transformation(
              extent={{-8,-8},{8,8}},
              rotation=0,
              origin={-40,-90})));

        parameter Modelica.SIunits.PressureDifference pump_dp_nominal = 1000000
          "Nominal pressure difference of pump at full speed";

        parameter Real pump_m_flow_max=5
          "Maximum mass flow generated by pump at full speed";

        parameter Real pump_speed_y_set = 0.8
          "Relative speed of pump (a value of 1 will generate maximum mass flow)";

      initial equation
        assert(pump_speed_y_set >= 0 and pump_speed_y_set <= 1,
         "Relative speed of pump must be between 0 and 1");

      equation
        connect(val_m_flow.port_a, senMasFlo.port_b)
          annotation (Line(points={{64,-70},{56,-70}}, color={0,127,255}));
        connect(m_flow_return_secondary_set, PID_m_flow.u_s)
          annotation (Line(points={{80,120},{80,57.2}}, color={0,0,127}));
        connect(jun_p_low.port_2, pump.port_a)
          annotation (Line(points={{-66,-70},{-26,-70}},color={0,127,255}));
        connect(jun_p_high.port_1, pump.port_b)
          annotation (Line(points={{6,-70},{-14,-70}},color={0,127,255}));
        connect(exp.port_a, pump.port_a) annotation (Line(points={{-40,-56},{-40,-70},
                {-26,-70}},color={0,127,255}));
        connect(jun_p_low.port_3, val_delta_p.port_a)
          annotation (Line(points={{-70,-66},{-70,-10},{-16,-10}},
                                                               color={0,127,255}));
        connect(val_delta_p.port_b, jun_p_high.port_3)
          annotation (Line(points={{-4,-10},{10,-10},{10,-66}},
                                                             color={0,127,255}));
        connect(port_a, pipe_return.port_a)
          annotation (Line(points={{-100,0},{-80,0},{-80,-14}}, color={0,127,255}));
        connect(port_b, pipe_supply.ports_b[1])
          annotation (Line(points={{100,0},{80,0},{80,-14}}, color={0,127,255}));
        connect(pipe_return.heatPort, ambientTemperature.port) annotation (Line(
              points={{-74,-20},{-50,-20},{-50,4}},                  color={191,0,0}));
        connect(pipe_supply.heatPort, ambientTemperature.port)
          annotation (Line(points={{74,-20},{-50,-20},{-50,4}},color={191,0,0}));
        connect(pipe_supply.port_a, val_m_flow.port_b) annotation (Line(points={{80,-26},
                {80,-70},{76,-70}},                             color={0,127,
                255}));
        connect(senRelPre.port_a, pump.port_b) annotation (Line(points={{-24,-30},{0,-30},
                {0,-70},{-14,-70}},color={0,127,255}));
        connect(senRelPre.port_b, pump.port_a) annotation (Line(points={{-36,-30},{-60,
                -30},{-60,-70},{-26,-70}},color={0,127,255}));
        connect(pump.m_flow_in, pump_speed_y.y)
          annotation (Line(points={{-20,-77.2},{-20,-90},{-31.2,-90}},
                                                                 color={0,0,127}));
        connect(jun_p_low.port_1, pipe_return.ports_b[1]) annotation (Line(points={{-74,
                -70},{-80,-70},{-80,-26}}, color={0,127,255}));
        connect(senMasFlo.port_a, coo.port_b)
          annotation (Line(points={{44,-70},{36,-70}}, color={0,127,255}));
        connect(coo.port_a, jun_p_high.port_2)
          annotation (Line(points={{24,-70},{14,-70}}, color={0,127,255}));
        connect(delta_p_secondary_set, PID_delta_p.u_s) annotation (Line(points=
               {{0,120},{0,57.2},{1.33227e-15,57.2}}, color={0,0,127}));
        connect(T_return_secondary_set, coo.TSet) annotation (Line(points={{-80,
                120},{-80,30},{20,30},{20,-65.2},{22.8,-65.2}}, color={0,0,127}));
        connect(PID_delta_p.y, val_delta_p.y) annotation (Line(points={{
                -1.11022e-15,43.4},{-1.11022e-15,10},{-10,10},{-10,-2.8}},
              color={0,0,127}));
        connect(senRelPre.p_rel, PID_delta_p.u_m) annotation (Line(points={{-30,
                -24.6},{-30,50},{-7.2,50}}, color={0,0,127}));
        connect(senMasFlo.m_flow, PID_m_flow.u_m) annotation (Line(points={{50,
                -63.4},{50,50},{72.8,50}}, color={0,0,127}));
        connect(PID_m_flow.y, val_m_flow.y) annotation (Line(points={{80,43.4},
                {80,10},{70,10},{70,-62.8}}, color={0,0,127}));
        annotation (Icon(graphics={Rectangle(
                extent={{-100,100},{100,-100}},
                lineColor={0,0,0},
                fillColor={51,100,0},
                fillPattern=FillPattern.Solid,
                radius=10)}));
      end SimplifiedPIDPressureCtrl;

    end SecondarySide;

    package BaseClasses
      partial model TeststandBase

        parameter Modelica.SIunits.MassFlowRate m1_flow_nominal "Nominal mass flow rate primary side"
          annotation(Dialog(group="Nominal condition"));
        parameter Modelica.SIunits.MassFlowRate m2_flow_nominal "Nominal mass flow rate secondary side"
          annotation(Dialog(group="Nominal condition"));
        parameter Modelica.SIunits.PressureDifference dp1_nominal(displayUnit="bar") "Nominal pressure difference primary side"
          annotation(Dialog(group="Nominal condition"));
        parameter Modelica.SIunits.PressureDifference dp2_nominal(displayUnit="bar") "Nominal pressure difference secondary side"
          annotation(Dialog(group="Nominal condition"));
        parameter Modelica.SIunits.PressureDifference dpValve_nominal "Nominal pressure drop of fully open valve"
          annotation(Dialog(group="Nominal condition"));

        replaceable BaseClasses.SubstationBase substation(
          m1_flow_nominal=m1_flow_nominal,
          m2_flow_nominal=m2_flow_nominal,
          dp1_nominal=dp1_nominal,
          dp2_nominal=dp2_nominal,
          dpValve_nominal=dpValve_nominal)
          annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
        Modelica.Blocks.Interfaces.RealInput delta_p_primary_set annotation (
            Placement(transformation(
              extent={{20,-20},{-20,20}},
              rotation=0,
              origin={120,60})));
        Modelica.Blocks.Interfaces.RealInput m_flow_return_secondary
          annotation (Placement(transformation(
              extent={{20,-20},{-20,20}},
              rotation=0,
              origin={120,0})));
        Modelica.Blocks.Interfaces.RealInput T_supply_primary_set annotation (
            Placement(transformation(
              extent={{-20,-20},{20,20}},
              rotation=0,
              origin={-120,60})));
        Modelica.Blocks.Interfaces.RealInput T_supply_secondary_set annotation (
           Placement(transformation(
              extent={{-20,-20},{20,20}},
              rotation=0,
              origin={-120,0})));
        replaceable PrimarySideBase primarySide(m_flow_nominal=m1_flow_nominal)
                                                annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=180,
              origin={0,60})));
        replaceable SecondarySideBase secondarySide(m_flow_nominal=m2_flow_nominal)
          annotation (Placement(transformation(extent={{-10,-50},{10,-70}})));
        Modelica.Blocks.Interfaces.RealInput T_return_secondary_set annotation (
           Placement(transformation(
              extent={{-20,-20},{20,20}},
              rotation=0,
              origin={-120,-60})));
        Modelica.Blocks.Interfaces.RealInput delta_p_secondary_set annotation (
            Placement(transformation(
              extent={{20,-20},{-20,20}},
              rotation=0,
              origin={120,-60})));
        IBPSA.Fluid.Sensors.TemperatureTwoPort sen_T_return_primary(redeclare
            package Medium = IBPSA.Media.Water, m_flow_nominal=primarySide.m_flow_nominal)
          annotation (Placement(transformation(
              extent={{-8,8},{8,-8}},
              rotation=90,
              origin={30,40})));
        IBPSA.Fluid.Sensors.MassFlowRate sen_m_flow_return_primary(redeclare
            package Medium = IBPSA.Media.Water) annotation (Placement(
              transformation(
              extent={{-8,-8},{8,8}},
              rotation=90,
              origin={30,20})));
        IBPSA.Fluid.Sensors.TemperatureTwoPort sen_T_supply_secondary(
            redeclare package Medium = IBPSA.Media.Water, m_flow_nominal=
              secondarySide.m_flow_nominal) annotation (Placement(
              transformation(
              extent={{8,8},{-8,-8}},
              rotation=90,
              origin={-30,-40})));
        IBPSA.Fluid.Sensors.MassFlowRate sen_m_flow_supply_primary(redeclare
            package Medium = IBPSA.Media.Water) annotation (Placement(
              transformation(
              extent={{8,-8},{-8,8}},
              rotation=90,
              origin={-30,-20})));
        Modelica.Blocks.Interfaces.RealOutput T_return_primary annotation (
            Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=90,
              origin={60,110})));
        Modelica.Blocks.Interfaces.RealOutput m_flow_return_primary annotation (
           Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=90,
              origin={-60,110})));
        Modelica.Blocks.Interfaces.RealOutput T_supply_secondary annotation (
            Placement(transformation(
              extent={{10,-10},{-10,10}},
              rotation=90,
              origin={60,-110})));
        Modelica.Blocks.Interfaces.RealOutput m_flow_supply_secondary annotation (
            Placement(transformation(
              extent={{10,-10},{-10,10}},
              rotation=90,
              origin={-60,-110})));
      equation
        connect(T_supply_primary_set, primarySide.T_supply_primary_set)
          annotation (Line(points={{-120,60},{-80,60},{-80,80},{-4,80},{-4,72}},
              color={0,0,127}));
        connect(delta_p_primary_set, primarySide.delta_p_primary_set)
          annotation (Line(points={{120,60},{80,60},{80,80},{4,80},{4,72}},
              color={0,0,127}));
        connect(T_supply_secondary_set, substation.T_supply_secondary_set)
          annotation (Line(points={{-120,0},{-12,0}}, color={0,0,127}));
        connect(delta_p_secondary_set, secondarySide.delta_p_secondary_set)
          annotation (Line(points={{120,-60},{80,-60},{80,-90},{0,-90},{0,-72}},
              color={0,0,127}));
        connect(T_return_secondary_set, secondarySide.T_return_secondary_set)
          annotation (Line(points={{-120,-60},{-80,-60},{-80,-90},{-8,-90},{-8,
                -72}}, color={0,0,127}));
        connect(secondarySide.m_flow_return_secondary_set,
          m_flow_return_secondary) annotation (Line(points={{8,-72},{8,-80},{70,
                -80},{70,0},{120,0}}, color={0,0,127}));
        connect(sen_m_flow_return_primary.port_b, sen_T_return_primary.port_a)
          annotation (Line(points={{30,28},{30,32}}, color={0,127,255}));
        connect(primarySide.port_a, sen_T_return_primary.port_b) annotation (
            Line(points={{10,60},{30,60},{30,48}}, color={0,127,255}));
        connect(substation.port_b1, sen_m_flow_return_primary.port_a)
          annotation (Line(points={{10,6},{30,6},{30,12}}, color={0,127,255}));
        connect(substation.port_a1, primarySide.port_b) annotation (Line(points=
               {{-10,6},{-30,6},{-30,60},{-10,60}}, color={0,127,255}));
        connect(sen_m_flow_return_primary.m_flow, m_flow_return_primary)
          annotation (Line(points={{21.2,20},{-60,20},{-60,110}}, color={0,0,
                127}));
        connect(sen_T_return_primary.T, T_return_primary) annotation (Line(
              points={{38.8,40},{60,40},{60,110}}, color={0,0,127}));
        connect(sen_m_flow_supply_primary.port_a, substation.port_b2)
          annotation (Line(points={{-30,-12},{-30,-6},{-10,-6}}, color={0,127,
                255}));
        connect(sen_m_flow_supply_primary.port_b, sen_T_supply_secondary.port_a)
          annotation (Line(points={{-30,-28},{-30,-32}}, color={0,127,255}));
        connect(secondarySide.port_a, sen_T_supply_secondary.port_b)
          annotation (Line(points={{-10,-60},{-30,-60},{-30,-48}}, color={0,127,
                255}));
        connect(substation.port_a2, secondarySide.port_b) annotation (Line(
              points={{10,-6},{30,-6},{30,-60},{10,-60}}, color={0,127,255}));
        connect(sen_m_flow_supply_primary.m_flow, m_flow_supply_secondary)
          annotation (Line(points={{-38.8,-20},{-60,-20},{-60,-110}}, color={0,0,127}));
        connect(sen_T_supply_secondary.T, T_supply_secondary) annotation (Line(
              points={{-21.2,-40},{60,-40},{60,-110}}, color={0,0,127}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
              Text(
                extent={{-149,-108},{151,-148}},
                lineColor={0,0,255},
                textString="%name")}),                                 Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end TeststandBase;

      partial model SubstationBase
        extends IBPSA.Fluid.Interfaces.PartialFourPortInterface(
          redeclare package Medium1 = IBPSA.Media.Water,
          redeclare package Medium2 = IBPSA.Media.Water);

        parameter Modelica.SIunits.PressureDifference dp1_nominal(displayUnit="bar") "Pressure difference" annotation(Dialog(group="Nominal condition"));
        parameter Modelica.SIunits.PressureDifference dp2_nominal(displayUnit="bar") "Pressure difference" annotation(Dialog(group="Nominal condition"));
        parameter Modelica.SIunits.PressureDifference dpValve_nominal(displayUnit="bar") "Nominal pressure drop of fully open valve" annotation(Dialog(group="Nominal condition"));

        Modelica.Blocks.Interfaces.RealInput T_supply_secondary_set
          annotation (Placement(transformation(extent={{-140,-20},{-100,20}})));
      end SubstationBase;

      partial model PrimarySideBase
        extends IBPSA.Fluid.Interfaces.PartialTwoPortInterface(
          redeclare package Medium = IBPSA.Media.Water);
        Modelica.Blocks.Interfaces.RealInput delta_p_primary_set annotation (
            Placement(transformation(extent={{-20,-20},{20,20}},
              rotation=90,
              origin={-40,-120})));
        Modelica.Blocks.Interfaces.RealInput T_supply_primary_set
          annotation (Placement(transformation(extent={{20,-20},{-20,20}},
              rotation=-90,
              origin={40,-120})));
      end PrimarySideBase;

      partial model SecondarySideBase
        extends IBPSA.Fluid.Interfaces.PartialTwoPortInterface(
          redeclare package Medium = IBPSA.Media.Water);
        Modelica.Blocks.Interfaces.RealInput m_flow_return_secondary_set
          annotation (Placement(transformation(
              extent={{-20,20},{20,-20}},
              rotation=-90,
              origin={80,120})));
        Modelica.Blocks.Interfaces.RealInput T_return_secondary_set annotation (
           Placement(transformation(
              extent={{-20,20},{20,-20}},
              rotation=270,
              origin={-80,120})));
        Modelica.Blocks.Interfaces.RealInput delta_p_secondary_set annotation (
            Placement(transformation(
              extent={{-20,20},{20,-20}},
              rotation=270,
              origin={0,120})));
      end SecondarySideBase;
    end BaseClasses;
  end Components;

  package Examples
    model TestStandDummy_Test
      TestStandDummy teststand(
        m1_flow_nominal=5,
        m2_flow_nominal=5,
        dp1_nominal=200000,
        dp2_nominal=200000,
        dpValve_nominal=200000,substation(show_T=true))
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
      Modelica.Blocks.Sources.Constant m_flow_secondary(k=3)
        annotation (Placement(transformation(extent={{80,-10},{60,10}})));
      Modelica.Blocks.Sources.Step T_secondary_set_degC(
        height=5,
        offset=40,
        startTime=100)
        annotation (Placement(transformation(extent={{-80,-60},{-60,-40}})));
      Modelica.Blocks.Sources.Constant T_supply_primary_set_degC(k=80)
        annotation (Placement(transformation(extent={{-80,40},{-60,60}})));
      Modelica.Blocks.Sources.Constant delta_p_primary_set_bar(k=6)
        annotation (Placement(transformation(extent={{80,40},{60,60}})));
      Modelica.Blocks.Sources.Constant T_supply_secondary_set_degC(k=60)
        annotation (Placement(transformation(extent={{-80,-10},{-60,10}})));
      Modelica.Blocks.Sources.Constant delta_p_secondary_set_bar(k=2)
        annotation (Placement(transformation(extent={{80,-60},{60,-40}})));
      Modelica.Blocks.Math.UnitConversions.From_degC T_supply_primary_set
        annotation (Placement(transformation(extent={{-50,40},{-30,60}})));
      Modelica.Blocks.Math.UnitConversions.From_bar delta_p_primary_set
        annotation (Placement(transformation(extent={{50,40},{30,60}})));
      Modelica.Blocks.Math.UnitConversions.From_degC T_supply_secondary_set
        annotation (Placement(transformation(extent={{-50,-10},{-30,10}})));
      Modelica.Blocks.Math.UnitConversions.From_degC T_secondary_set
        annotation (Placement(transformation(extent={{-50,-60},{-30,-40}})));
      Modelica.Blocks.Math.UnitConversions.From_bar delta_p_secondary_set
        annotation (Placement(transformation(extent={{50,-60},{30,-40}})));
    equation
      connect(teststand.m_flow_return_secondary, m_flow_secondary.y)
        annotation (Line(points={{12,0},{59,0}}, color={0,0,127}));
      connect(teststand.T_supply_primary_set, T_supply_primary_set.y)
        annotation (Line(points={{-12,6},{-20,6},{-20,50},{-29,50}}, color={0,0,
              127}));
      connect(T_supply_primary_set.u, T_supply_primary_set_degC.y)
        annotation (Line(points={{-52,50},{-59,50}}, color={0,0,127}));
      connect(teststand.T_supply_secondary_set, T_supply_secondary_set.y)
        annotation (Line(points={{-12,0},{-29,0}}, color={0,0,127}));
      connect(T_supply_secondary_set.u, T_supply_secondary_set_degC.y)
        annotation (Line(points={{-52,0},{-59,0}}, color={0,0,127}));
      connect(teststand.T_return_secondary_set, T_secondary_set.y) annotation (
          Line(points={{-12,-6},{-20,-6},{-20,-50},{-29,-50}}, color={0,0,127}));
      connect(T_secondary_set.u, T_secondary_set_degC.y)
        annotation (Line(points={{-52,-50},{-59,-50}}, color={0,0,127}));
      connect(teststand.delta_p_primary_set, delta_p_primary_set.y) annotation (
         Line(points={{12,6},{20,6},{20,50},{29,50}}, color={0,0,127}));
      connect(delta_p_primary_set.u, delta_p_primary_set_bar.y)
        annotation (Line(points={{52,50},{59,50}}, color={0,0,127}));
      connect(teststand.delta_p_secondary_set, delta_p_secondary_set.y)
        annotation (Line(points={{12,-6},{20,-6},{20,-50},{29,-50}}, color={0,0,
              127}));
      connect(delta_p_secondary_set.u, delta_p_secondary_set_bar.y)
        annotation (Line(points={{52,-50},{59,-50}}, color={0,0,127}));
      annotation (
        Icon(coordinateSystem(preserveAspectRatio=false), graphics={
            Ellipse(lineColor = {75,138,73},
                    fillColor={255,255,255},
                    fillPattern = FillPattern.Solid,
                    extent={{-100,-100},{100,100}}),
            Polygon(lineColor = {0,0,255},
                    fillColor = {75,138,73},
                    pattern = LinePattern.None,
                    fillPattern = FillPattern.Solid,
                    points={{-36,60},{64,0},{-36,-60},{-36,60}})}),
        Diagram(coordinateSystem(preserveAspectRatio=false)),
        experiment(StopTime=300, __Dymola_Algorithm="Dassl"));
    end TestStandDummy_Test;

    model TestStandSimplifiedDirectPressureCtrl_Test
      TestStandSimplifiedDirectPressureCtrl
                                         teststand(
        m1_flow_nominal=5,
        m2_flow_nominal=5,
        dp1_nominal=200000,
        dp2_nominal=200000,
        dpValve_nominal=200000)
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
      Modelica.Blocks.Sources.Constant m_flow_secondary(k=3)
        annotation (Placement(transformation(extent={{80,-10},{60,10}})));
      Modelica.Blocks.Sources.Step T_return_secondary_set_degC(
        height=5,
        offset=40,
        startTime=100)
        annotation (Placement(transformation(extent={{-80,-60},{-60,-40}})));
      Modelica.Blocks.Sources.Constant T_supply_primary_set_degC(k=80)
        annotation (Placement(transformation(extent={{-80,40},{-60,60}})));
      Modelica.Blocks.Sources.Step     delta_p_primary_set_bar(
        height=-0.5,
        offset=6,
        startTime=200)
        annotation (Placement(transformation(extent={{80,40},{60,60}})));
      Modelica.Blocks.Sources.Constant T_supply_secondary_set_degC(k=60)
        annotation (Placement(transformation(extent={{-80,-10},{-60,10}})));
      Modelica.Blocks.Sources.Step     delta_p_secondary_set_bar(
        height=-0.4,
        offset=2,
        startTime=150)
        annotation (Placement(transformation(extent={{80,-60},{60,-40}})));
      Modelica.Blocks.Math.UnitConversions.From_degC T_supply_primary_set
        annotation (Placement(transformation(extent={{-50,40},{-30,60}})));
      Modelica.Blocks.Math.UnitConversions.From_degC T_supply_secondary_set
        annotation (Placement(transformation(extent={{-50,-10},{-30,10}})));
      Modelica.Blocks.Math.UnitConversions.From_degC T_secondary_set
        annotation (Placement(transformation(extent={{-50,-60},{-30,-40}})));
      Modelica.Blocks.Math.UnitConversions.From_bar delta_p_primary_set
        annotation (Placement(transformation(extent={{50,40},{30,60}})));
      Modelica.Blocks.Math.UnitConversions.From_bar delta_p_secondary_set
        annotation (Placement(transformation(extent={{50,-60},{30,-40}})));
    equation
      connect(teststand.m_flow_return_secondary, m_flow_secondary.y)
        annotation (Line(points={{12,0},{59,0}}, color={0,0,127}));
      connect(teststand.delta_p_primary_set, delta_p_primary_set.y) annotation (
         Line(points={{12,6},{20,6},{20,50},{29,50}}, color={0,0,127}));
      connect(delta_p_primary_set.u, delta_p_primary_set_bar.y)
        annotation (Line(points={{52,50},{59,50}}, color={0,0,127}));
      connect(teststand.delta_p_secondary_set, delta_p_secondary_set.y)
        annotation (Line(points={{12,-6},{20,-6},{20,-50},{29,-50}}, color={0,0,
              127}));
      connect(delta_p_secondary_set.u, delta_p_secondary_set_bar.y)
        annotation (Line(points={{52,-50},{59,-50}}, color={0,0,127}));
      connect(teststand.T_return_secondary_set, T_secondary_set.y) annotation (
          Line(points={{-12,-6},{-20,-6},{-20,-50},{-29,-50}}, color={0,0,127}));
      connect(T_secondary_set.u, T_return_secondary_set_degC.y)
        annotation (Line(points={{-52,-50},{-59,-50}}, color={0,0,127}));
      connect(teststand.T_supply_secondary_set, T_supply_secondary_set.y)
        annotation (Line(points={{-12,0},{-29,0}}, color={0,0,127}));
      connect(T_supply_secondary_set.u, T_supply_secondary_set_degC.y)
        annotation (Line(points={{-52,0},{-59,0}}, color={0,0,127}));
      connect(teststand.T_supply_primary_set, T_supply_primary_set.y)
        annotation (Line(points={{-12,6},{-20,6},{-20,50},{-29,50}}, color={0,0,
              127}));
      connect(T_supply_primary_set.u, T_supply_primary_set_degC.y)
        annotation (Line(points={{-52,50},{-59,50}}, color={0,0,127}));
      annotation (
        Icon(coordinateSystem(preserveAspectRatio=false), graphics={
            Ellipse(lineColor = {75,138,73},
                    fillColor={255,255,255},
                    fillPattern = FillPattern.Solid,
                    extent={{-100,-100},{100,100}}),
            Polygon(lineColor = {0,0,255},
                    fillColor = {75,138,73},
                    pattern = LinePattern.None,
                    fillPattern = FillPattern.Solid,
                    points={{-36,60},{64,0},{-36,-60},{-36,60}})}),
        Diagram(coordinateSystem(preserveAspectRatio=false)),
        experiment(StopTime=300, __Dymola_Algorithm="Dassl"));
    end TestStandSimplifiedDirectPressureCtrl_Test;

    model TestStandSimplifiedPIDPressureCtrl_Test
      TestStandSimplifiedPIDPressureCtrl teststand(
        m1_flow_nominal=5,
        m2_flow_nominal=5,
        dp1_nominal=200000,
        dp2_nominal=200000,
        dpValve_nominal=200000)
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
      Modelica.Blocks.Sources.Constant m_flow_secondary(k=3)
        annotation (Placement(transformation(extent={{80,-10},{60,10}})));
      Modelica.Blocks.Sources.Step T_return_secondary_set_degC(
        height=5,
        offset=40,
        startTime=100)
        annotation (Placement(transformation(extent={{-80,-60},{-60,-40}})));
      Modelica.Blocks.Sources.Constant T_supply_primary_set_degC(k=80)
        annotation (Placement(transformation(extent={{-80,40},{-60,60}})));
      Modelica.Blocks.Sources.Step     delta_p_primary_set_bar(
        height=-0.5,
        offset=6,
        startTime=200)
        annotation (Placement(transformation(extent={{80,40},{60,60}})));
      Modelica.Blocks.Sources.Constant T_supply_secondary_set_degC(k=60)
        annotation (Placement(transformation(extent={{-80,-10},{-60,10}})));
      Modelica.Blocks.Sources.Step     delta_p_secondary_set_bar(
        height=-0.4,
        offset=2,
        startTime=150)
        annotation (Placement(transformation(extent={{80,-60},{60,-40}})));
      Modelica.Blocks.Math.UnitConversions.From_degC T_supply_primary_set
        annotation (Placement(transformation(extent={{-50,40},{-30,60}})));
      Modelica.Blocks.Math.UnitConversions.From_degC T_supply_secondary_set
        annotation (Placement(transformation(extent={{-50,-10},{-30,10}})));
      Modelica.Blocks.Math.UnitConversions.From_degC T_secondary_set
        annotation (Placement(transformation(extent={{-50,-60},{-30,-40}})));
      Modelica.Blocks.Math.UnitConversions.From_bar delta_p_primary_set
        annotation (Placement(transformation(extent={{50,40},{30,60}})));
      Modelica.Blocks.Math.UnitConversions.From_bar delta_p_secondary_set
        annotation (Placement(transformation(extent={{50,-60},{30,-40}})));
    equation
      connect(teststand.m_flow_return_secondary, m_flow_secondary.y)
        annotation (Line(points={{12,0},{59,0}}, color={0,0,127}));
      connect(teststand.T_supply_primary_set, T_supply_primary_set.y)
        annotation (Line(points={{-12,6},{-20,6},{-20,50},{-29,50}}, color={0,0,
              127}));
      connect(T_supply_primary_set.u, T_supply_primary_set_degC.y)
        annotation (Line(points={{-52,50},{-59,50}}, color={0,0,127}));
      connect(teststand.T_supply_secondary_set, T_supply_secondary_set.y)
        annotation (Line(points={{-12,0},{-29,0}}, color={0,0,127}));
      connect(T_supply_secondary_set.u, T_supply_secondary_set_degC.y)
        annotation (Line(points={{-52,0},{-59,0}}, color={0,0,127}));
      connect(teststand.T_return_secondary_set, T_secondary_set.y) annotation (
          Line(points={{-12,-6},{-20,-6},{-20,-50},{-29,-50}}, color={0,0,127}));
      connect(T_secondary_set.u, T_return_secondary_set_degC.y)
        annotation (Line(points={{-52,-50},{-59,-50}}, color={0,0,127}));
      connect(teststand.delta_p_primary_set, delta_p_primary_set.y) annotation (
         Line(points={{12,6},{20,6},{20,50},{29,50}}, color={0,0,127}));
      connect(delta_p_primary_set.u, delta_p_primary_set_bar.y)
        annotation (Line(points={{52,50},{59,50}}, color={0,0,127}));
      connect(teststand.delta_p_secondary_set, delta_p_secondary_set.y)
        annotation (Line(points={{12,-6},{20,-6},{20,-50},{29,-50}}, color={0,0,
              127}));
      connect(delta_p_secondary_set.u, delta_p_secondary_set_bar.y)
        annotation (Line(points={{52,-50},{59,-50}}, color={0,0,127}));
      annotation (
        Icon(coordinateSystem(preserveAspectRatio=false), graphics={
            Ellipse(lineColor = {75,138,73},
                    fillColor={255,255,255},
                    fillPattern = FillPattern.Solid,
                    extent={{-100,-100},{100,100}}),
            Polygon(lineColor = {0,0,255},
                    fillColor = {75,138,73},
                    pattern = LinePattern.None,
                    fillPattern = FillPattern.Solid,
                    points={{-36,60},{64,0},{-36,-60},{-36,60}})}),
        Diagram(coordinateSystem(preserveAspectRatio=false)),
        experiment(StopTime=300, __Dymola_Algorithm="Dassl"));
    end TestStandSimplifiedPIDPressureCtrl_Test;
    annotation (Icon(graphics={
          Rectangle(
            lineColor={200,200,200},
            fillColor={248,248,248},
            fillPattern=FillPattern.HorizontalCylinder,
            extent={{-100,-100},{100,100}},
            radius=25),
          Polygon(
            origin={20,14},
            lineColor={78,138,73},
            fillColor={78,138,73},
            pattern=LinePattern.None,
            fillPattern=FillPattern.Solid,
            points={{-58.0,46.0},{42.0,-14.0},{-58.0,-74.0},{-58.0,46.0}}),
          Rectangle(
            extent={{-100,100},{100,-100}},
            lineColor={135,135,135},
            radius=25)}));
  end Examples;

  package ModelCalibration
    model CalibrateHeatExchanger

      IBPSA.Fluid.Sources.MassFlowSource_T source_secondary(
        redeclare package Medium = IBPSA.Media.Water,
        use_m_flow_in=true,
        use_T_in=true,
        nPorts=1)
        annotation (Placement(transformation(extent={{-50,-30},{-30,-10}})));
      IBPSA.Fluid.Sources.MassFlowSource_T source_primary(
        use_T_in=true,
        redeclare package Medium = IBPSA.Media.Water,
        use_m_flow_in=true,
        nPorts=1) annotation (Placement(transformation(extent={{50,10},{30,30}})));
      IBPSA.Fluid.Sources.Boundary_pT sink_primary(
        redeclare package Medium = IBPSA.Media.Water,
        use_T_in=false,
        nPorts=1) annotation (Placement(transformation(extent={{-50,10},{-30,30}})));
      IBPSA.Fluid.Sources.Boundary_pT sink_secondary(
        redeclare package Medium = IBPSA.Media.Water,
        use_T_in=false,
        nPorts=1) annotation (Placement(transformation(extent={{50,-30},{30,-10}})));
      IBPSA.Fluid.HeatExchangers.ConstantEffectiveness hex(
        redeclare package Medium1 = IBPSA.Media.Water,
        redeclare package Medium2 = IBPSA.Media.Water,
        m1_flow_nominal=m1_flow_nominal,
        m2_flow_nominal=m2_flow_nominal,
        dp1_nominal=dp1_nominal,
        dp2_nominal=dp2_nominal,
        eps=eps)
        annotation (Placement(transformation(extent={{10,-10},{-10,10}})));
      IBPSA.Fluid.Sensors.Temperature T_primary_return(redeclare package Medium =
            IBPSA.Media.Water)
        annotation (Placement(transformation(extent={{-30,40},{-10,60}})));
      IBPSA.Fluid.Sensors.Temperature T_secondary_supply(redeclare package
          Medium =
            IBPSA.Media.Water)
        annotation (Placement(transformation(extent={{30,-40},{10,-60}})));
      Modelica.Blocks.Interfaces.RealInput m_flow_secondary_return
        annotation (Placement(transformation(extent={{-140,20},{-100,60}})));
      Modelica.Blocks.Interfaces.RealInput T_secondary_return_degC
        annotation (Placement(transformation(extent={{-140,-60},{-100,-20}})));
      Modelica.Blocks.Math.UnitConversions.From_degC T_secondary_return_K
        annotation (Placement(transformation(extent={{-90,-50},{-70,-30}})));
      Modelica.Blocks.Interfaces.RealInput m_flow_primary_supply
        annotation (Placement(transformation(extent={{140,20},{100,60}})));
      Modelica.Blocks.Interfaces.RealInput T_primary_supply_degC
        annotation (Placement(transformation(extent={{140,-60},{100,-20}})));
      Modelica.Blocks.Math.UnitConversions.From_degC T_primary_supply_K
        annotation (Placement(transformation(extent={{90,-50},{70,-30}})));
      parameter Modelica.SIunits.Efficiency eps=0.9 "Heat exchanger effectiveness";
      parameter Modelica.SIunits.MassFlowRate m1_flow_nominal=4
        "Nominal mass flow rate";
      parameter Modelica.SIunits.MassFlowRate m2_flow_nominal=4
        "Nominal mass flow rate";
      parameter Modelica.SIunits.PressureDifference dp1_nominal=2000
        "Pressure difference";
      parameter Modelica.SIunits.PressureDifference dp2_nominal=500
        "Pressure difference";
    equation
      connect(hex.port_a1, source_primary.ports[1]) annotation (Line(points={{10,
              6},{20,6},{20,20},{30,20}}, color={0,127,255}));
      connect(hex.port_b1, sink_primary.ports[1]) annotation (Line(points={{-10,6},
              {-20,6},{-20,20},{-30,20}}, color={0,127,255}));
      connect(hex.port_b2, sink_secondary.ports[1]) annotation (Line(points={{10,
              -6},{20,-6},{20,-20},{30,-20}}, color={0,127,255}));
      connect(hex.port_a2, source_secondary.ports[1]) annotation (Line(points={{
              -10,-6},{-20,-6},{-20,-20},{-30,-20}}, color={0,127,255}));
      connect(hex.port_b1, T_primary_return.port)
        annotation (Line(points={{-10,6},{-20,6},{-20,40}}, color={0,127,255}));
      connect(hex.port_b2, T_secondary_supply.port)
        annotation (Line(points={{10,-6},{20,-6},{20,-40}}, color={0,127,255}));
      connect(m_flow_secondary_return, source_secondary.m_flow_in) annotation (Line(
            points={{-120,40},{-60,40},{-60,-12},{-52,-12}}, color={0,0,127}));
      connect(T_secondary_return_degC, T_secondary_return_K.u)
        annotation (Line(points={{-120,-40},{-92,-40}}, color={0,0,127}));
      connect(T_secondary_return_K.y, source_secondary.T_in) annotation (Line(
            points={{-69,-40},{-60,-40},{-60,-16},{-52,-16}}, color={0,0,127}));
      connect(m_flow_primary_supply, source_primary.m_flow_in) annotation (Line(
            points={{120,40},{60,40},{60,28},{52,28}}, color={0,0,127}));
      connect(T_primary_supply_degC, T_primary_supply_K.u)
        annotation (Line(points={{120,-40},{92,-40}}, color={0,0,127}));
      connect(T_primary_supply_K.y, source_primary.T_in) annotation (Line(points={{69,
              -40},{60,-40},{60,24},{52,24}}, color={0,0,127}));
      annotation (                    experiment(StopTime=300, __Dymola_Algorithm=
              "Dassl"), Icon(graphics={
            Rectangle(
              lineColor={28,108,200},
              fillColor={248,248,248},
              fillPattern=FillPattern.HorizontalCylinder,
              extent={{-100,-100},{100,100}},
              radius=10),
            Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={0,0,0},
              radius=10)}));
    end CalibrateHeatExchanger;

    model CalibrateSubstationPIDMassFlow

      IBPSA.Fluid.Sources.MassFlowSource_T source_secondary(
        redeclare package Medium = IBPSA.Media.Water,
        use_m_flow_in=true,
        use_T_in=true,
        nPorts=1)
        annotation (Placement(transformation(extent={{-60,-30},{-40,-10}})));
      IBPSA.Fluid.Sources.Boundary_pT sink_primary(
        redeclare package Medium = IBPSA.Media.Water,
        use_T_in=false,
        T=333.15,
        nPorts=1) annotation (Placement(transformation(extent={{-60,10},{-40,30}})));
      IBPSA.Fluid.Sources.Boundary_pT sink_secondary(
        redeclare package Medium = IBPSA.Media.Water,
        use_T_in=false,
        nPorts=1) annotation (Placement(transformation(extent={{62,-80},{42,-60}})));
      IBPSA.Fluid.HeatExchangers.ConstantEffectiveness hex(
        redeclare package Medium1 = IBPSA.Media.Water,
        redeclare package Medium2 = IBPSA.Media.Water,
        m1_flow_nominal=4,
        m2_flow_nominal=4,
        dp1_nominal(displayUnit="bar") = 400000,
        dp2_nominal(displayUnit="bar") = 100000,
        eps=0.9)
        annotation (Placement(transformation(extent={{10,-10},{-10,10}})));
      Modelica.Blocks.Math.UnitConversions.From_degC T_return_secondary_K
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=90,
            origin={-100,-90})));
      IBPSA.Fluid.Sources.Boundary_pT source_primary(
        redeclare package Medium = IBPSA.Media.Water,
        use_p_in=true,
        use_T_in=true,
        nPorts=1)
        annotation (Placement(transformation(extent={{180,40},{160,60}})));
      IBPSA.Fluid.Actuators.Valves.TwoWayLinear val(
        redeclare package Medium = IBPSA.Media.Water,
        m_flow_nominal=4,
        dpValve_nominal(displayUnit="bar") = 100000,
        riseTime=10,
        y_start=0,
        dpFixed_nominal(displayUnit="bar") = 50000)
        annotation (Placement(transformation(extent={{62,60},{42,40}})));
      Modelica.Blocks.Continuous.LimPID PID(
        controllerType=Modelica.Blocks.Types.SimpleController.PID,
        k=k,
        Ti=Ti,
        Td=Td,
        yMax=1,
        yMin=0) annotation (Placement(transformation(extent={{100,10},{80,30}})));
      IBPSA.Fluid.Sensors.Temperature T_supply_secondary_sense_K(redeclare
          package Medium = IBPSA.Media.Water)
        annotation (Placement(transformation(extent={{40,-20},{60,-40}})));
      Modelica.Blocks.Math.UnitConversions.From_degC T_supply_secondary_set_K
        annotation (Placement(transformation(extent={{160,-10},{140,10}})));
      Modelica.Blocks.Math.UnitConversions.From_bar delta_p_supply_primary_Pa
        annotation (Placement(transformation(extent={{120,90},{140,110}})));
      Modelica.Blocks.Sources.Constant p_supply_primary_baseline_Pa(k=IBPSA.Media.Water.p_default)
        annotation (Placement(transformation(extent={{80,70},{100,90}})));
      Modelica.Blocks.Math.Add add
        annotation (Placement(transformation(extent={{160,100},{180,80}})));
      Modelica.Blocks.Math.UnitConversions.From_degC T_supply_primary_K
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=90,
            origin={100,-90})));
      parameter Real k=0.005 "Gain of controller";
      parameter Modelica.SIunits.Time Ti=0.5
        "Time constant of Integrator block";
      parameter Modelica.SIunits.Time Td=0.1
        "Time constant of Derivative block";
      Modelica.Blocks.Interfaces.RealInput delta_p_supply_primary_bar
        annotation (Placement(transformation(
            extent={{-20,-20},{20,20}},
            rotation=-90,
            origin={100,140})));
      Modelica.Blocks.Interfaces.RealInput T_supply_primary_degC annotation (
          Placement(transformation(
            extent={{-20,-20},{20,20}},
            rotation=90,
            origin={100,-140})));
      Modelica.Blocks.Interfaces.RealInput T_return_secondary_degC annotation (
          Placement(transformation(
            extent={{-20,-20},{20,20}},
            rotation=90,
            origin={-100,-140})));
      Modelica.Blocks.Interfaces.RealInput m_flow_supply_secondary annotation (
          Placement(transformation(
            extent={{-20,-20},{20,20}},
            rotation=270,
            origin={-100,140})));
      Modelica.Blocks.Interfaces.RealInput T_supply_secondary_set_degC
        annotation (Placement(transformation(
            extent={{-20,20},{20,-20}},
            rotation=180,
            origin={220,0})));
    equation
      connect(hex.port_b1, sink_primary.ports[1]) annotation (Line(points={{-10,6},
              {-20,6},{-20,20},{-40,20}}, color={0,127,255}));
      connect(hex.port_b2, sink_secondary.ports[1]) annotation (Line(points={{10,-6},
              {20,-6},{20,-70},{42,-70}},     color={0,127,255}));
      connect(hex.port_a2, source_secondary.ports[1]) annotation (Line(points={{-10,-6},
              {-20,-6},{-20,-20},{-40,-20}},         color={0,127,255}));
      connect(source_secondary.T_in, T_return_secondary_K.y) annotation (Line(
            points={{-62,-16},{-100,-16},{-100,-79}}, color={0,0,127}));
      connect(hex.port_a1, val.port_b) annotation (Line(points={{10,6},{20,6},{
              20,50},{42,50}},
                             color={0,127,255}));
      connect(val.port_a, source_primary.ports[1])
        annotation (Line(points={{62,50},{160,50}},  color={0,127,255}));
      connect(PID.y, val.y)
        annotation (Line(points={{79,20},{52,20},{52,38}},    color={0,0,127}));
      connect(T_supply_secondary_sense_K.port, hex.port_b2)
        annotation (Line(points={{50,-20},{50,-6},{10,-6}}, color={0,127,255}));
      connect(T_supply_secondary_set_K.y, PID.u_s)
        annotation (Line(points={{139,0},{120,0},{120,20},{102,20}},
                                                   color={0,0,127}));
      connect(T_supply_secondary_sense_K.T, PID.u_m)
        annotation (Line(points={{57,-30},{90,-30},{90,8}},   color={0,0,127}));
      connect(p_supply_primary_baseline_Pa.y, add.u1) annotation (Line(points={{101,80},
              {150,80},{150,84},{158,84}},           color={0,0,127}));
      connect(delta_p_supply_primary_Pa.y, add.u2) annotation (Line(points={{
              141,100},{152,100},{152,96},{158,96}}, color={0,0,127}));
      connect(add.y, source_primary.p_in) annotation (Line(points={{181,90},{190,
              90},{190,58},{182,58}}, color={0,0,127}));
      connect(T_supply_primary_K.y, source_primary.T_in) annotation (Line(
            points={{100,-79},{100,-60},{190,-60},{190,54},{182,54}}, color={0,
              0,127}));
      connect(delta_p_supply_primary_Pa.u, delta_p_supply_primary_bar)
        annotation (Line(points={{118,100},{100,100},{100,140}}, color={0,0,127}));
      connect(T_supply_primary_degC, T_supply_primary_K.u)
        annotation (Line(points={{100,-140},{100,-102}}, color={0,0,127}));
      connect(T_return_secondary_degC, T_return_secondary_K.u)
        annotation (Line(points={{-100,-140},{-100,-102}}, color={0,0,127}));
      connect(m_flow_supply_secondary, source_secondary.m_flow_in) annotation (
          Line(points={{-100,140},{-100,-12},{-62,-12}}, color={0,0,127}));
      connect(T_supply_secondary_set_K.u, T_supply_secondary_set_degC)
        annotation (Line(points={{162,0},{220,0}}, color={0,0,127}));
      annotation (                    experiment(StopTime=300, __Dymola_Algorithm=
              "Dassl"),
        Diagram(coordinateSystem(extent={{-200,-120},{200,120}})),
        Icon(coordinateSystem(extent={{-200,-120},{200,120}}), graphics={
            Rectangle(
              lineColor={28,108,200},
              fillColor={248,248,248},
              fillPattern=FillPattern.HorizontalCylinder,
              extent={{-200,-120},{200,120}},
              radius=10),
            Rectangle(
              extent={{-200,120},{200,-120}},
              lineColor={0,0,0},
              radius=10)}));
    end CalibrateSubstationPIDMassFlow;

    model CalibrateSecondarySidePIDMassFlow

      IBPSA.Fluid.Sources.Boundary_pT sink_primary(
        redeclare package Medium = IBPSA.Media.Water,
        use_T_in=false,
        T=313.15,
        nPorts=1) annotation (Placement(transformation(extent={{-60,10},{-40,30}})));
      IBPSA.Fluid.Sources.Boundary_pT sink_secondary(
        redeclare package Medium = IBPSA.Media.Water,
        use_T_in=false,
        nPorts=1) annotation (Placement(transformation(extent={{62,-80},{42,-60}})));
      Modelica.Blocks.Continuous.LimPID PID_primary(
        controllerType=Modelica.Blocks.Types.SimpleController.PID,
        k=0.005,
        Ti=0.5,
        Td=0.1,
        yMax=1,
        yMin=0)
        annotation (Placement(transformation(extent={{100,10},{120,30}})));
      IBPSA.Fluid.HeatExchangers.ConstantEffectiveness hex(
        redeclare package Medium1 = IBPSA.Media.Water,
        redeclare package Medium2 = IBPSA.Media.Water,
        m1_flow_nominal=4,
        m2_flow_nominal=4,
        dp1_nominal(displayUnit="bar") = 400000,
        dp2_nominal(displayUnit="bar") = 100000,
        eps=0.9)
        annotation (Placement(transformation(extent={{10,-10},{-10,10}})));
      Modelica.Blocks.Math.UnitConversions.From_degC T_return_secondary_K
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=90,
            origin={-100,-90})));
      IBPSA.Fluid.Sources.Boundary_pT source_primary(
        redeclare package Medium = IBPSA.Media.Water,
        use_p_in=true,
        use_T_in=true,
        nPorts=1)
        annotation (Placement(transformation(extent={{180,40},{160,60}})));
      IBPSA.Fluid.Actuators.Valves.TwoWayLinear val(
        redeclare package Medium = IBPSA.Media.Water,
        m_flow_nominal=4,
        dpValve_nominal(displayUnit="bar") = 100000,
        riseTime=10,
        y_start=0,
        dpFixed_nominal(displayUnit="bar") = 50000)
        annotation (Placement(transformation(extent={{140,60},{120,40}})));
      IBPSA.Fluid.Sensors.Temperature T_supply_secondary_sense_K(redeclare
          package Medium = IBPSA.Media.Water)
        annotation (Placement(transformation(extent={{40,-20},{60,-40}})));
      Modelica.Blocks.Sources.Constant T_supply_secondary_set_degC(k=50)
        annotation (Placement(transformation(extent={{40,10},{60,30}})));
      Modelica.Blocks.Math.UnitConversions.From_degC T_supply_secondary_set_K
        annotation (Placement(transformation(extent={{70,10},{90,30}})));
      Modelica.Blocks.Math.UnitConversions.From_bar delta_p_supply_primary_Pa
        annotation (Placement(transformation(extent={{120,90},{140,110}})));
      Modelica.Blocks.Sources.Constant p_supply_primary_baseline_Pa(k=IBPSA.Media.Water.p_default)
        annotation (Placement(transformation(extent={{80,70},{100,90}})));
      Modelica.Blocks.Math.Add p_supply_primary_Pa
        annotation (Placement(transformation(extent={{160,100},{180,80}})));
      Modelica.Blocks.Math.UnitConversions.From_degC T_supply_primary_K annotation (
         Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=90,
            origin={100,-90})));
      Modelica.Blocks.Interfaces.RealInput delta_p_supply_primary_bar annotation (
          Placement(transformation(
            extent={{-20,-20},{20,20}},
            rotation=-90,
            origin={100,140})));
      Modelica.Blocks.Interfaces.RealInput T_supply_primary_degC annotation (
          Placement(transformation(
            extent={{-20,-20},{20,20}},
            rotation=90,
            origin={100,-140})));
      Modelica.Blocks.Interfaces.RealInput T_return_secondary_degC annotation (
          Placement(transformation(
            extent={{-20,-20},{20,20}},
            rotation=90,
            origin={-100,-140})));
      Modelica.Blocks.Interfaces.RealInput delta_p_return_secondary_bar annotation (
         Placement(transformation(
            extent={{-20,-20},{20,20}},
            rotation=270,
            origin={-100,140})));
      Modelica.Blocks.Sources.Constant m_flow_secondary(k=3)
        annotation (Placement(transformation(extent={{-130,-30},{-110,-10}})));
      IBPSA.Fluid.Sources.Boundary_pT source_secondary1(
        redeclare package Medium = IBPSA.Media.Water,
        use_p_in=true,
        use_T_in=true,
        nPorts=1)
        annotation (Placement(transformation(extent={{-160,-60},{-140,-40}})));
      Modelica.Blocks.Math.UnitConversions.From_bar delta_p_return_secondary_Pa
        annotation (Placement(transformation(extent={{-120,90},{-140,110}})));
      Modelica.Blocks.Sources.Constant p_return_secondary_baseline_Pa(k=IBPSA.Media.Water.p_default)
        annotation (Placement(transformation(extent={{-80,70},{-100,90}})));
      Modelica.Blocks.Math.Add p_return_secondary_bar
        annotation (Placement(transformation(extent={{-160,100},{-180,80}})));
      IBPSA.Fluid.Sensors.MassFlowRate senMasFlo(redeclare package Medium =
            IBPSA.Media.Water)
        annotation (Placement(transformation(extent={{-100,-60},{-80,-40}})));
      IBPSA.Fluid.Actuators.Valves.TwoWayLinear val_secondary(
        redeclare package Medium = IBPSA.Media.Water,
        m_flow_nominal=4,
        dpValve_nominal(displayUnit="bar") = 100000,
        riseTime=10,
        y_start=0,
        dpFixed_nominal(displayUnit="bar") = 50000)
        annotation (Placement(transformation(extent={{-60,-60},{-40,-40}})));
      Modelica.Blocks.Continuous.LimPID PID_secondary(
        controllerType=Modelica.Blocks.Types.SimpleController.PID,
        k=k,
        Ti=Ti,
        Td=Td,
        yMax=1,
        yMin=0)
        annotation (Placement(transformation(extent={{-100,-30},{-80,-10}})));
      parameter Real k=0.03 "Gain of controller";
      parameter Modelica.SIunits.Time Ti=0.5
        "Time constant of Integrator block";
      parameter Modelica.SIunits.Time Td=0.1
        "Time constant of Derivative block";
    equation
      connect(hex.port_b1, sink_primary.ports[1]) annotation (Line(points={{-10,6},
              {-20,6},{-20,20},{-40,20}}, color={0,127,255}));
      connect(hex.port_b2, sink_secondary.ports[1]) annotation (Line(points={{10,-6},
              {20,-6},{20,-70},{42,-70}},     color={0,127,255}));
      connect(hex.port_a1, val.port_b) annotation (Line(points={{10,6},{20,6},{20,
              50},{120,50}}, color={0,127,255}));
      connect(val.port_a, source_primary.ports[1])
        annotation (Line(points={{140,50},{160,50}}, color={0,127,255}));
      connect(PID_primary.y, val.y) annotation (Line(points={{121,20},{130,20},
              {130,38}}, color={0,0,127}));
      connect(T_supply_secondary_sense_K.port, hex.port_b2)
        annotation (Line(points={{50,-20},{50,-6},{10,-6}}, color={0,127,255}));
      connect(T_supply_secondary_set_degC.y, T_supply_secondary_set_K.u)
        annotation (Line(points={{61,20},{68,20}}, color={0,0,127}));
      connect(T_supply_secondary_set_K.y, PID_primary.u_s)
        annotation (Line(points={{91,20},{98,20}}, color={0,0,127}));
      connect(T_supply_secondary_sense_K.T, PID_primary.u_m) annotation (Line(
            points={{57,-30},{110,-30},{110,8}}, color={0,0,127}));
      connect(p_supply_primary_baseline_Pa.y, p_supply_primary_Pa.u1) annotation (
          Line(points={{101,80},{150,80},{150,84},{158,84}}, color={0,0,127}));
      connect(delta_p_supply_primary_Pa.y, p_supply_primary_Pa.u2) annotation (Line(
            points={{141,100},{150,100},{150,96},{158,96}}, color={0,0,127}));
      connect(p_supply_primary_Pa.y, source_primary.p_in) annotation (Line(points={{
              181,90},{190,90},{190,58},{182,58}}, color={0,0,127}));
      connect(T_supply_primary_K.y, source_primary.T_in) annotation (Line(points={{100,
              -79},{100,-60},{190,-60},{190,54},{182,54}}, color={0,0,127}));
      connect(delta_p_supply_primary_Pa.u, delta_p_supply_primary_bar)
        annotation (Line(points={{118,100},{100,100},{100,140}}, color={0,0,127}));
      connect(T_supply_primary_degC, T_supply_primary_K.u)
        annotation (Line(points={{100,-140},{100,-102}}, color={0,0,127}));
      connect(T_return_secondary_degC, T_return_secondary_K.u)
        annotation (Line(points={{-100,-140},{-100,-102}}, color={0,0,127}));
      connect(senMasFlo.m_flow,PID_secondary. u_m)
        annotation (Line(points={{-90,-39},{-90,-32}}, color={0,0,127}));
      connect(PID_secondary.y,val_secondary. y) annotation (Line(points={{-79,-20},{
              -50,-20},{-50,-38}},  color={0,0,127}));
      connect(PID_secondary.u_s,m_flow_secondary. y)
        annotation (Line(points={{-102,-20},{-109,-20}},
                                                       color={0,0,127}));
      connect(source_secondary1.ports[1], senMasFlo.port_a)
        annotation (Line(points={{-140,-50},{-100,-50}}, color={0,127,255}));
      connect(val_secondary.port_a,senMasFlo. port_b)
        annotation (Line(points={{-60,-50},{-80,-50}}, color={0,127,255}));
      connect(val_secondary.port_b, hex.port_a2) annotation (Line(points={{-40,-50},
              {-20,-50},{-20,-6},{-10,-6}},      color={0,127,255}));
      connect(p_return_secondary_bar.u2, delta_p_return_secondary_Pa.y) annotation (
         Line(points={{-158,96},{-148,96},{-148,100},{-141,100}}, color={0,0,127}));
      connect(p_return_secondary_bar.y, source_secondary1.p_in) annotation (Line(
            points={{-181,90},{-190,90},{-190,-42},{-162,-42}}, color={0,0,127}));
      connect(p_return_secondary_bar.u1,p_return_secondary_baseline_Pa. y)
        annotation (Line(points={{-158,84},{-148,84},{-148,80},{-101,80}},
            color={0,0,127}));
      connect(T_return_secondary_K.y, source_secondary1.T_in) annotation (Line(
            points={{-100,-79},{-100,-70},{-190,-70},{-190,-46},{-162,-46}}, color={
              0,0,127}));
      connect(delta_p_return_secondary_Pa.u, delta_p_return_secondary_bar)
        annotation (Line(points={{-118,100},{-100,100},{-100,140}}, color={0,0,127}));
      annotation (                    experiment(StopTime=300, __Dymola_Algorithm=
              "Dassl"),
        Diagram(coordinateSystem(extent={{-200,-120},{200,120}})),
        Icon(coordinateSystem(extent={{-200,-120},{200,120}}), graphics={
            Rectangle(
              lineColor={28,108,200},
              fillColor={248,248,248},
              fillPattern=FillPattern.HorizontalCylinder,
              extent={{-200,-120},{200,120}},
              radius=10),
            Rectangle(
              extent={{-200,120},{200,-120}},
              lineColor={0,0,0},
              radius=10)}));
    end CalibrateSecondarySidePIDMassFlow;
  end ModelCalibration;
  annotation (uses(IBPSA(version="3.0.0"), Modelica(version="3.2.3")));
end SubstationTeststand;
