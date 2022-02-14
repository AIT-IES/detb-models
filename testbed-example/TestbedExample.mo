within ;
package TestbedExample

  package Components
    "This sub-package contains components used for modelling the thermal system of the multi-energy network benchmark"

    model AggregatedConsumers
      "Aggregated consumer in the thermal network, separated from the supply line via substation."

      IBPSA.Fluid.Sensors.RelativePressure senRelPre(redeclare package Medium =
            IBPSA.Media.Water)
        annotation (Placement(transformation(extent={{-10,-70},{10,-50}})));
      extends IBPSA.Fluid.Interfaces.PartialTwoPort(
        redeclare package Medium = IBPSA.Media.Water);

      parameter Modelica.SIunits.MassFlowRate m_flow_nominal = 10 "Nominal mass flow rate"
           annotation (Dialog(group="Nominal conditions"));
      parameter Modelica.SIunits.PressureDifference dp_nominal(displayUnit="bar")=
           1000000 "Nominal pressure difference"
               annotation (Dialog(group="Nominal conditions"));
      parameter String fileName="modelica://TestbedExample/data/network_data.txt" "File where matrix is stored"
           annotation(Dialog(group="Aggregated demand profiles"),loadSelector(filter="Text files(*.txt);;CSV files (*.csv)",caption="Open data file"));
      parameter String tableName="NetworkData" "Table name on file"
           annotation (Dialog(group="Aggregated demand profiles"));
      parameter Integer columns[:]={2,3,4} "Columns of table to be interpolated"
          annotation (Dialog(group="Aggregated demand profiles"));

      Modelica.Blocks.Interfaces.RealOutput dp annotation (Placement(
            transformation(
            extent={{10,-10},{-10,10}},
            rotation=90,
            origin={0,-110})));
      Modelica.Blocks.Sources.CombiTimeTable profiles(
        tableOnFile=true,
        tableName=tableName,
        columns=columns,
        extrapolation=Modelica.Blocks.Types.Extrapolation.HoldLastPoint,
        fileName=Modelica.Utilities.Files.loadResource(fileName))
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=90,
            origin={0,10})));
      DisHeatLib.Utilities.ThermoHydraulicEquivalent.DemandWithDelay_mQ
        demand_mQ(
        m_flow_nominal=m_flow_nominal,
        dp_nominal=dp_nominal,
        T_min=328.15,
        pipeType=pipeType,
        length=length,
        fac=fac,
        cf=cf) annotation (Placement(transformation(extent={{-10,70},{10,50}})));
      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a port_ht
        annotation (Placement(transformation(extent={{-10,90},{10,110}})));

      replaceable parameter DisHeatLib.Pipes.BaseClasses.BasePipe pipeType
        "Pipe type"
         annotation (choicesAllMatching=true, Dialog(group="Plug Flow Delay"));
      parameter Modelica.SIunits.Length length=1000 "Effective pipe length"
        annotation (Dialog(group="Plug Flow Delay"));
      parameter Real fac=1
        "Factor to take into account flow resistance of bends etc., fac=dp_nominal/dpStraightPipe_nominal"
        annotation (Dialog(group="Plug Flow Delay"));
      parameter Real cf=7
        "Correction factor of heat losses (needed for aggregation)"
        annotation (Dialog(group="Plug Flow Delay"));
    equation
      connect(demand_mQ.m_flow_set,profiles. y[2])
        annotation (Line(points={{-4,48},{-4,30},{-8.88178e-16,30},{-8.88178e-16,21}},
                                  color={0,0,127}));
      connect(demand_mQ.Q_demand_set,profiles. y[1])
        annotation (Line(points={{4,48},{4,30},{-8.88178e-16,30},{-8.88178e-16,21}},
                                  color={0,0,127}));
      connect(demand_mQ.port_a, port_a) annotation (Line(points={{-10,60},{-60,60},{
              -60,0},{-100,0}}, color={0,127,255}));
      connect(demand_mQ.port_b, port_b) annotation (Line(points={{10,60},{60,60},{60,
              0},{100,0}}, color={0,127,255}));
      connect(senRelPre.port_b, port_b) annotation (Line(points={{10,-60},{60,-60},{
              60,0},{100,0}}, color={0,127,255}));
      connect(senRelPre.port_a, port_a) annotation (Line(points={{-10,-60},{-60,-60},
              {-60,0},{-100,0}}, color={0,127,255}));
      connect(senRelPre.p_rel, dp)
        annotation (Line(points={{0,-69},{0,-110}}, color={0,0,127}));
      connect(port_ht, demand_mQ.port_ht)
        annotation (Line(points={{0,100},{0,86},{0,70},{6,70}},
                                                  color={191,0,0}));
      annotation (
        Icon(
          coordinateSystem(preserveAspectRatio=false),
          graphics={
            Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={0,0,0},
              fillColor={227,221,16},
              fillPattern=FillPattern.Solid,
              radius=15),
            Ellipse(
              extent={{-88,88},{88,-88}},
              lineColor={0,0,0},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),
            Bitmap(
              extent={{-64,-46},{68,32}}, fileName="modelica://TestbedExample/images/consumer.svg"),
            Polygon(
              points={{-132,92},{-132,92}},
              lineColor={0,0,0},
              lineThickness=0.5,
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid)}),
        Diagram(coordinateSystem(preserveAspectRatio=false)));
    end AggregatedConsumers;

    model SingleConsumer
      "Single consumer in the thermal network, connected via a substation."

      extends IBPSA.Fluid.Interfaces.PartialTwoPort(redeclare package Medium =
            IBPSA.Media.Water);

      parameter Boolean use_p_ref=true "If true, include an expansion vessel that acts as an internal reference point for the absolute pressure level."
        annotation(Evaluate=true, HideResult=true, Dialog(tab="Advanced"));

      parameter Modelica.SIunits.MassFlowRate m_flow_nominal "Nominal mass flow rate secondary side"
        annotation(Dialog(group="Nominal condition"));
      parameter Modelica.SIunits.PressureDifference dp_nominal(displayUnit=
            "bar") "Nominal pressure difference secondary side"
            annotation(Dialog(group="Nominal condition"));

      parameter Modelica.SIunits.Temperature T_supply_nominal "Nominal secondary supply temperature" annotation(Dialog(group="Nominal condition"));
      parameter Modelica.SIunits.Temperature T_return_nominal "Nominal secondary return temperature" annotation(Dialog(group="Nominal condition"));

      parameter Modelica.SIunits.Power Q_flow_nominal(displayUnit="kW") = 10000
        "Nominal heat flow rate" annotation(Dialog(group="Nominal condition"));

      parameter String fileName="modelica://TestbedExample/data/space_heating_load_profile_single_building_1week.txt"
        "File where matrix is stored"
            annotation(Dialog(group="Demand profile"),loadSelector(filter="Text files(*.txt);;CSV files (*.csv)",caption="Open data file"));
      parameter String tableName="HeatDemand"
        "Table name on file"
            annotation(Dialog(group="Demand profile"));
      parameter Real scaling=2 "Scaling factor for heat demand"
        annotation(Dialog(group="Demand profile"));


      DisHeatLib.Demand.Demand demand(
        redeclare package Medium = IBPSA.Media.Water,
        allowFlowReversal=true,
        show_T=true,
        from_dp=true,
        dp_nominal(displayUnit="bar") = dp_nominal,
        Q_flow_nominal(displayUnit="kW") = Q_flow_nominal,
        TemSup_nominal=T_supply_nominal,
        TemRet_nominal=T_return_nominal,
        heatLoad=DisHeatLib.Demand.BaseClasses.InputTypeDemand.FileQ,
        scaling=scaling,
        tableName=tableName,
        fileName=fileName,
        redeclare DisHeatLib.Demand.BaseDemands.Radiator demandType)
        annotation (Placement(transformation(extent={{-10,-10},{10,10}},
            rotation=0,
            origin={50,0})));
      IBPSA.Fluid.Storage.ExpansionVessel exp(redeclare package Medium =
            IBPSA.Media.Water, V_start=2) if use_p_ref
        annotation (Placement(transformation(extent={{-60,20},{-40,40}})));

      IBPSA.Fluid.Movers.FlowControlled_dp fan(
        redeclare package Medium = IBPSA.Media.Water,
        m_flow_nominal=m_flow_nominal,
        redeclare IBPSA.Fluid.Movers.Data.Pumps.Wilo.Stratos25slash1to4 per,
        inputType=IBPSA.Fluid.Types.InputType.Constant,
        addPowerToMedium=false,
        dp_nominal=dp_nominal)   annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={0,0})));
    equation

      connect(fan.port_b, demand.port_a)
        annotation (Line(points={{10,0},{40,0}},           color={0,127,255}));
      connect(exp.port_a, fan.port_a) annotation (Line(points={{-50,20},{-50,0},{-10,
              0}},                         color={0,127,255}));
      connect(demand.port_b, port_b)
        annotation (Line(points={{60,0},{100,0}}, color={0,127,255}));
      connect(port_a, fan.port_a)
        annotation (Line(points={{-100,0},{-10,0}}, color={0,127,255}));
      annotation (
        Icon(
          coordinateSystem(preserveAspectRatio=false),
          graphics={
            Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={0,0,0},
              fillColor={227,221,16},
              fillPattern=FillPattern.Solid,
              radius=15),
            Ellipse(
              extent={{-88,88},{88,-88}},
              lineColor={0,0,0},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),
            Bitmap(
              extent={{-54,-56},{54,52}}, fileName="modelica://TestbedExample/images/heat-exchanger.svg"),
            Polygon(
              points={{-132,92},{-132,92}},
              lineColor={0,0,0},
              lineThickness=0.5,
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid)}),
        Diagram(coordinateSystem(preserveAspectRatio=false)));
    end SingleConsumer;

    model SingleConsumerWithBooster
      "Single consumer in the thermal network, connected via a substation."

      extends IBPSA.Fluid.Interfaces.PartialTwoPort(redeclare package Medium =
            IBPSA.Media.Water);

      parameter Boolean use_p_ref=true "If true, include an expansion vessel that acts as an internal reference point for the absolute pressure level."
        annotation(Evaluate=true, HideResult=true, Dialog(tab="Advanced"));

      parameter Modelica.SIunits.MassFlowRate m_flow_nominal "Nominal mass flow rate secondary side"
        annotation(Dialog(group="Nominal condition"));
      parameter Modelica.SIunits.PressureDifference dp_nominal(displayUnit=
            "bar") "Nominal pressure difference secondary side"
            annotation(Dialog(group="Nominal condition"));

      parameter Modelica.SIunits.Temperature T_supply_nominal "Nominal secondary supply temperature" annotation(Dialog(group="Nominal condition"));
      parameter Modelica.SIunits.Temperature T_return_nominal "Nominal secondary return temperature" annotation(Dialog(group="Nominal condition"));

      parameter Modelica.SIunits.Power Q_flow_nominal(displayUnit="kW") = 10000
        "Nominal heat flow rate" annotation(Dialog(group="Nominal condition"));

      parameter String fileName="modelica://TestbedExample/data/space_heating_load_profile_single_building_1week.txt"
        "File where matrix is stored"
            annotation(Dialog(group="Demand profile"),loadSelector(filter="Text files(*.txt);;CSV files (*.csv)",caption="Open data file"));
      parameter String tableName="HeatDemand"
        "Table name on file"
            annotation(Dialog(group="Demand profile"));
      parameter Real scaling=2 "Scaling factor for heat demand"
        annotation(Dialog(group="Demand profile"));

      DisHeatLib.Demand.Demand demand(
        redeclare package Medium = IBPSA.Media.Water,
        allowFlowReversal=true,
        show_T=true,
        from_dp=true,
        dp_nominal(displayUnit="bar") = dp_nominal,
        Q_flow_nominal(displayUnit="kW") = Q_flow_nominal,
        TemSup_nominal=T_supply_nominal,
        TemRet_nominal=T_return_nominal,
        heatLoad=DisHeatLib.Demand.BaseClasses.InputTypeDemand.FileQ,
        scaling=scaling,
        tableName=tableName,
        fileName=fileName,
        redeclare DisHeatLib.Demand.BaseDemands.Radiator demandType)
        annotation (Placement(transformation(extent={{-10,-10},{10,10}},
            rotation=0,
            origin={62,0})));
      IBPSA.Fluid.Storage.ExpansionVessel exp(redeclare package Medium =
            IBPSA.Media.Water, V_start=2) if use_p_ref
        annotation (Placement(transformation(extent={{-80,20},{-60,40}})));

      IBPSA.Fluid.Movers.FlowControlled_dp fan(
        redeclare package Medium = IBPSA.Media.Water,
        m_flow_nominal=m_flow_nominal,
        redeclare IBPSA.Fluid.Movers.Data.Pumps.Wilo.Stratos25slash1to4 per,
        inputType=IBPSA.Fluid.Types.InputType.Constant,
        addPowerToMedium=false,
        dp_nominal=dp_nominal)   annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={20,0})));
      Heatpump heatpump(
        m1_flow_nominal=m_flow_nominal,
        m2_flow_nominal=m_flow_nominal,
        Q_flow_evaporator_rated=Q_flow_nominal,
        P_rated(displayUnit="kW") = Q_flow_nominal,
        T_cond_out_max=363.15,
        T_cond_out_target=T_supply_nominal,
        evaporator(m_flow_nominal=0.05, dp_nominal=100000))
        annotation (Placement(transformation(extent={{-40,-30},{-20,-10}})));
      IBPSA.Fluid.Sources.FixedBoundary hpEvapSink(
        redeclare package Medium = IBPSA.Media.Water,
        p=300000,
        T=293.15,
        nPorts=1)
        annotation (Placement(transformation(extent={{-80,-50},{-60,-30}})));
      IBPSA.Fluid.Sources.FixedBoundary hpEvapSource(
        redeclare package Medium = IBPSA.Media.Water,
        p=400000,
        T=313.15,
        nPorts=1)
        annotation (Placement(transformation(extent={{20,-50},{0,-30}})));
    equation

      connect(fan.port_b, demand.port_a)
        annotation (Line(points={{30,0},{52,0}},           color={0,127,255}));
      connect(demand.port_b, port_b)
        annotation (Line(points={{72,0},{100,0}}, color={0,127,255}));
      connect(port_a, heatpump.port_a1) annotation (Line(points={{-100,0},{-50,
              0},{-50,-14},{-40,-14}}, color={0,127,255}));
      connect(exp.port_a, heatpump.port_a1) annotation (Line(points={{-70,20},{
              -70,0},{-50,0},{-50,-14},{-40,-14}}, color={0,127,255}));
      connect(heatpump.port_b1, fan.port_a) annotation (Line(points={{-20,-14},
              {-10,-14},{-10,0},{10,0}}, color={0,127,255}));
      connect(hpEvapSink.ports[1], heatpump.port_b2) annotation (Line(points={{
              -60,-40},{-50,-40},{-50,-26},{-40,-26}}, color={0,127,255}));
      connect(hpEvapSource.ports[1], heatpump.port_a2) annotation (Line(points=
              {{0,-40},{-10,-40},{-10,-26},{-20,-26}}, color={0,127,255}));
      annotation (
        Icon(
          coordinateSystem(preserveAspectRatio=false),
          graphics={
            Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={0,0,0},
              fillColor={244,125,35},
              fillPattern=FillPattern.Solid,
              radius=15),
            Ellipse(
              extent={{-88,88},{88,-88}},
              lineColor={0,0,0},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),
            Bitmap(
              extent={{-54,-56},{54,52}}, fileName="modelica://TestbedExample/images/heat-exchanger.svg"),
            Polygon(
              points={{-132,92},{-132,92}},
              lineColor={0,0,0},
              lineThickness=0.5,
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid)}),
        Diagram(coordinateSystem(preserveAspectRatio=false)));
    end SingleConsumerWithBooster;

    model ExternalGridConstantPressure
      "External thermal grid with constant pressure."
      extends IBPSA.Fluid.Interfaces.PartialTwoPortInterface(
        redeclare package Medium = IBPSA.Media.Water,
        m_flow_nominal=10);

      DisHeatLib.Supply.Supply_pT supply_pT(
        redeclare package Medium = IBPSA.Media.Water,
        Q_flow_nominal(displayUnit="MW") = Q_flow_nominal,
        TemSup_nominal=T_supply_nominal,
        TemRet_nominal=T_return_nominal,
        isElectric=false,
        powerCha(Q_flow={0}, P={0}),
        dp_nominal(displayUnit="bar") = dp_nominal,
        SupplyTemperature=DisHeatLib.Supply.BaseClasses.InputTypeSupplyTemp.Constant,
        TemOut_min=283.15,
        TemOut_max=283.15,
        TemSup_min=283.15,
        TemSup_max=283.15,
        dp_controller=false,
        dp_min=50000,
        dp_set=80000,
        dp_max=200000,
        ports_b,
        heater(show_T=true),
        nPorts=1)
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
      IBPSA.Fluid.Sensors.EntropyFlowRate sflow_supply(redeclare package Medium =
            IBPSA.Media.Water, m_flow_nominal=m_flow_nominal)
        annotation (Placement(transformation(extent={{60,-10},{80,10}})));
      IBPSA.Fluid.Sensors.Temperature T_supply(redeclare package Medium =
            IBPSA.Media.Water)
        annotation (Placement(transformation(extent={{30,20},{50,40}})));

      parameter Modelica.SIunits.Temperature T_supply_nominal(displayUnit="degC") = 358.15 "Nominal supply temperature";
      parameter Modelica.SIunits.Temperature T_return_nominal(displayUnit="degC") = 313.15 "Nominal return temperature";
      parameter Modelica.SIunits.HeatFlowRate Q_flow_nominal(displayUnit="kW") = 2000000 "Nominal heat flow rate";

      Modelica.SIunits.Heat Q_total_supply(start=0,displayUnit="kWh");
      Modelica.SIunits.Heat Q_total_return(start=0,displayUnit="kWh");

      IBPSA.Fluid.Sensors.EntropyFlowRate sflow_return(redeclare package Medium =
            IBPSA.Media.Water, m_flow_nominal=m_flow_nominal)
        annotation (Placement(transformation(extent={{-82,-10},{-62,10}})));
      IBPSA.Fluid.Sensors.Temperature T_return(redeclare package Medium =
            IBPSA.Media.Water)
        annotation (Placement(transformation(extent={{-50,20},{-30,40}})));
      parameter Modelica.SIunits.PressureDifference dp_nominal(displayUnit=
            "bar") = 600000 "Nominal pressure difference of pump";
    equation

      der(Q_total_supply) = T_supply.T*sflow_supply.S_flow;
      der(Q_total_return) = T_return.T*sflow_return.S_flow;

      connect(port_b,sflow_supply. port_b)
        annotation (Line(points={{100,0},{80,0}}, color={0,127,255}));
      connect(sflow_supply.port_a, supply_pT.ports_b[1])
        annotation (Line(points={{60,0},{10,0}},  color={0,127,255}));
      connect(sflow_supply.port_a, T_supply.port)
        annotation (Line(points={{60,0},{40,0},{40,20}},color={0,127,255}));
      connect(port_a, sflow_return.port_a)
        annotation (Line(points={{-100,0},{-82,0}}, color={0,127,255}));
      connect(sflow_return.port_b, supply_pT.port_a)
        annotation (Line(points={{-62,0},{-10,0}}, color={0,127,255}));
      connect(sflow_return.port_b, T_return.port)
        annotation (Line(points={{-62,0},{-40,0},{-40,20}}, color={0,127,255}));
      annotation (Icon(graphics={
            Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={0,0,0},
              fillColor={227,221,16},
              fillPattern=FillPattern.Solid,
              radius=15),
            Ellipse(
              extent={{-88,89},{88,-87}},
              lineColor={0,0,0},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),
            Bitmap(extent={{-56,-66},{58,8}}, fileName=
                  "modelica://TestbedExample/images/external-grid.svg"),
            Text(
              extent={{-62,88},{64,2}},
              lineColor={0,0,0},
              pattern=LinePattern.None,
              fillColor={0,0,0},
              fillPattern=FillPattern.None,
              textString="P")}));
    end ExternalGridConstantPressure;

    model Heatpump
      "Heat pump used in the power-to-heat facility."
      extends IBPSA.Fluid.Interfaces.PartialFourPortInterface(
      redeclare package Medium1 = IBPSA.Media.Water,
      redeclare package Medium2 = IBPSA.Media.Water);

      IBPSA.Fluid.HeatExchangers.HeaterCooler_u evaporator(
        redeclare package Medium = IBPSA.Media.Water,
        m_flow_nominal=10,
        show_T=false,
        dp_nominal=1000,
        T_start=333.15,
        Q_flow_nominal(displayUnit="kW") = -Q_flow_evaporator_rated)
        annotation (Placement(transformation(extent={{10,-40},{-10,-20}})));
      IBPSA.Fluid.HeatExchangers.Heater_T condenser(
        redeclare package Medium = IBPSA.Media.Water,
        m_flow_nominal=10,
        show_T=false,
        dp_nominal(displayUnit="bar") = 10000)
        annotation (Placement(transformation(extent={{-10,20},{10,40}})));
      IBPSA.Fluid.Sensors.Temperature TcondIn(redeclare package Medium =
            IBPSA.Media.Water)
        annotation (Placement(transformation(extent={{-86,46},{-74,34}})));
      IBPSA.Fluid.Sensors.Temperature TcondOut(redeclare package Medium =
            IBPSA.Media.Water)
        annotation (Placement(transformation(extent={{86,46},{74,34}})));
      IBPSA.Fluid.Sensors.Temperature TevapIn(redeclare package Medium =
            IBPSA.Media.Water)
        annotation (Placement(transformation(extent={{86,-46},{74,-34}})));
      IBPSA.Fluid.Sensors.Temperature TevapOut(redeclare package Medium =
            IBPSA.Media.Water)
        annotation (Placement(transformation(extent={{-86,-46},{-74,-34}})));
      Modelica.Blocks.Sources.RealExpression P_effective(y=P_0 + W_effective/eta_comp)
        annotation (Placement(transformation(extent={{-60,-90},{-20,-70}})));
      Modelica.Blocks.Interfaces.RealOutput P_el_heatpump annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=-90,
            origin={0,-110})));

      parameter Modelica.SIunits.HeatFlowRate Q_flow_evaporator_rated(displayUnit="kW") = 100000
        "Heat flow rating at evaporator";

      parameter Modelica.SIunits.Power P_rated(displayUnit="kW") = 100000.0 "Electrical power rating";
      parameter Modelica.SIunits.Power P_0(displayUnit="kW") = 300 "Electrical stand-by power consumption";

      parameter Modelica.SIunits.Temperature T_evap_out_min(displayUnit="degC") = 288.15 "Minimal evaporator outlet temperature";
      parameter Modelica.SIunits.Temperature T_cond_out_max(displayUnit="degC") = 358.15 "Maximum condenser outlet temperature";
      parameter Modelica.SIunits.Temperature T_cond_out_target(displayUnit="degC") = 348.15 "Condenser outlet temperature setpoint";

      parameter Modelica.SIunits.Efficiency eta_sys = 0.63 "Ratio between work provided by the pump and available thermodynamic work";
      parameter Modelica.SIunits.Efficiency eta_comp = 0.65 "Compressor efficiency";
      parameter Modelica.SIunits.DampingCoefficient lambda_comp = 0.05 "Compressor time constant";

      Modelica.SIunits.Power W_cond_max "Maximum mechanical work done in the condenser";
      Modelica.SIunits.Power W_evap_max "Maximum mechanical work done in the evaporator";
      Modelica.SIunits.Power W_max "Maximum mechanical work constraint";
      Modelica.SIunits.Power W_requested "Requested mechanical work";
      Modelica.SIunits.Power W_effective "Effective mechanical work";

      Modelica.SIunits.Temperature T_cond_log "Logarithmic mean temperature of condenser";
      Modelica.SIunits.Temperature T_evap_log "Logarithmic mean temperature of evaporator";
      Modelica.SIunits.Efficiency eta "Overall efficiency";

    protected
      parameter Modelica.SIunits.Power W_rated(fixed = false) "Mechanical power rating";
      parameter Modelica.SIunits.Power W_0(fixed = false) "Mechanical stand-by power consumption";
      parameter Modelica.SIunits.SpecificHeatCapacity cp_cond(fixed = false)
        "Specific heat capacity of medium 1 at nominal condition";
      parameter Modelica.SIunits.SpecificHeatCapacity cp_evap(fixed = false)
        "Specific heat capacity of medium 2 at nominal condition";

    initial equation
      W_rated = P_rated * eta_comp;
      W_0 = P_0 * eta_comp;

      cp_cond = Medium1.specificHeatCapacityCp(state_a1_inflow);
      cp_evap = Medium2.specificHeatCapacityCp(state_a2_inflow);

    equation
      T_cond_log = Utility.logMean(TcondIn.T, TcondOut.T);
      T_evap_log = Utility.logMean(TevapIn.T, TevapOut.T);

      eta = Utility.eta(T_evap_log, T_cond_log, eta_sys);

      W_cond_max = (T_cond_out_max - TcondIn.T) * cp_cond * m1_flow / eta;
      W_evap_max = (TevapIn.T - T_evap_out_min) * cp_evap * m2_flow / (eta - 1);
      W_max = max(0.0, min([W_evap_max, W_cond_max, W_rated]));
      W_requested = Utility.clamp(0, condenser.Q_flow/eta, W_max);

      evaporator.u = Utility.clamp(0, abs((1 - eta)*W_effective/Q_flow_evaporator_rated), W_effective/Q_flow_evaporator_rated);
      condenser.TSet = T_cond_out_target;

      der(W_effective) = lambda_comp * (W_requested - W_effective);

      connect(port_a2, evaporator.port_a)
        annotation (Line(points={{100,-60},{40,-60},{40,-30},{10,-30}},color={0,127,255}));
      connect(evaporator.port_b, port_b2)
        annotation (Line(points={{-10,-30},{-40,-30},{-40,-60},{-100,-60}},color={0,127,255}));
      connect(port_a1, condenser.port_a)
        annotation (Line(points={{-100,60},{-40,60},{-40,30},{-10,30}},color={0,127,255}));
      connect(port_b1, condenser.port_b)
        annotation (Line(points={{100,60},{40,60},{40,30},{10,30}},color={0,127,255}));
      connect(port_a1, TcondIn.port)
        annotation (Line(points={{-100,60},{-80,60},{-80,46}}, color={0,127,255}));
      connect(port_b1, TcondOut.port)
        annotation (Line(points={{100,60},{80,60},{80,46}}, color={0,127,255}));
      connect(TevapIn.port, port_a2)
        annotation (Line(points={{80,-46},{80,-60},{100,-60}}, color={0,127,255}));
      connect(port_b2, TevapOut.port) annotation (Line(points={{-100,-60},{-80,-60},
              {-80,-46}}, color={0,127,255}));

      connect(P_el_heatpump, P_effective.y)
        annotation (Line(points={{0,-110},{0,-80},{-18,-80}}, color={0,0,127}));

      annotation (Icon(graphics={
            Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={0,0,0},
              fillColor={227,221,16},
              fillPattern=FillPattern.Solid,
              radius=15),
            Ellipse(
              extent={{-88,88},{88,-88}},
              lineColor={0,0,0},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),
            Bitmap(extent={{-58,-60},{64,58}}, fileName=
                  "modelica://TestbedExample/images/heat-pump.svg"),
            Rectangle(
              extent={{-100,66},{0,54}},
              lineThickness=1,
              fillColor={28,108,200},
              fillPattern=FillPattern.Solid,
              pattern=LinePattern.None,
              lineColor={0,0,0}),
            Rectangle(
              extent={{0,-54},{100,-66}},
              lineThickness=1,
              fillColor={28,108,200},
              fillPattern=FillPattern.Solid,
              pattern=LinePattern.None,
              lineColor={0,0,0}),
            Rectangle(
              extent={{0,66},{100,54}},
              lineThickness=1,
              fillColor={238,46,47},
              fillPattern=FillPattern.Solid,
              pattern=LinePattern.None,
              lineColor={0,0,0}),
            Rectangle(
              extent={{-100,-54},{0,-66}},
              lineThickness=1,
              fillColor={0,0,173},
              fillPattern=FillPattern.Solid,
              pattern=LinePattern.None,
              lineColor={0,0,0})}));
    end Heatpump;

    model DHNetwork
      DisHeatLib.Boundary.SoilTemperature soilTemperature(
        inputType=DisHeatLib.Boundary.BaseClasses.InputTypeSoilTemp.Constant,
        T_const(displayUnit="degC") = 283.15,
        T_mean=283.15,
        T_amp(displayUnit="degC") = 10,
        t_min=0) annotation (Placement(transformation(extent={{40,60},{60,40}})));
      DisHeatLib.Pipes.DualPipe pipe1(
        show_T=true,
        redeclare package Medium = IBPSA.Media.Water,
        redeclare DisHeatLib.Pipes.Library.Isoplus.Isoplus_Std_DN100 pipeType,
        L=500,
        nPorts1=1,
        nPorts2=1)
        annotation (Placement(transformation(extent={{-40,-40},{-20,-20}})));
      DisHeatLib.Pipes.DualPipe pipe2(
        show_T=true,
        redeclare package Medium = IBPSA.Media.Water,
        redeclare DisHeatLib.Pipes.Library.Isoplus.Isoplus_Std_DN100 pipeType,
        L=10,
        nPorts2=1,
        nPorts1=1) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=90,
            origin={20,10})));
      IBPSA.Fluid.FixedResistances.Junction junction2(
        redeclare package Medium = IBPSA.Media.Water,
        m_flow_nominal={10,-20,10},
        dp_nominal={10,-10,10})
        annotation (Placement(transformation(extent={{40,-26},{20,-46}})));
      IBPSA.Fluid.FixedResistances.Junction junction1(
        redeclare package Medium = IBPSA.Media.Water,
        m_flow_nominal={20,-10,-10},
        dp_nominal={10,-10,-10})
        annotation (Placement(transformation(extent={{0,-14},{20,-34}})));
      Components.ExternalGridConstantPressure externalGrid(show_T=true)
        annotation (Placement(transformation(extent={{-80,-20},{-60,0}})));
      Components.AggregatedConsumers aggregatedConsumer(fileName=fileName,
        tableName=tableName,                            redeclare
          DisHeatLib.Pipes.Library.Isoplus.Isoplus_Std_DN50 pipeType)
        annotation (Placement(transformation(extent={{60,-20},{80,0}})));
      Modelica.Fluid.Interfaces.FluidPort_a port_a(redeclare package Medium =
            IBPSA.Media.Water)
        annotation (Placement(transformation(extent={{40,90},{60,110}})));
      Modelica.Fluid.Interfaces.FluidPort_b port_b(redeclare package Medium =
            IBPSA.Media.Water)
        annotation (Placement(transformation(extent={{-60,90},{-40,110}})));
      parameter String fileName="modelica://TestbedExample/data/network_data.txt"
        "File where matrix is stored";

      parameter String tableName="NetworkData" "Table name on file";

    equation
      connect(soilTemperature.port,pipe1. port_ht) annotation (Line(points={{50,40},
              {50,30},{-30,30},{-30,-20}}, color={191,0,0}));
      connect(junction1.port_1,pipe1. ports_b1[1])
        annotation (Line(points={{0,-24},{-20,-24}}, color={0,127,255}));
      connect(junction1.port_3,pipe2. port_a1) annotation (Line(points={{10,-14},
              {10,-8},{14,-8},{14,0}}, color={0,127,255}));
      connect(pipe2.port_ht,soilTemperature. port) annotation (Line(points={{10,10},
              {-30,10},{-30,30},{50,30},{50,40}},  color={191,0,0}));
      connect(pipe1.port_a2,junction2. port_2)
        annotation (Line(points={{-20,-36},{20,-36}}, color={0,127,255}));
      connect(pipe2.ports_b2[1],junction2. port_3) annotation (Line(points={{26,0},{
              26,-8},{30,-8},{30,-26}},    color={0,127,255}));
      connect(externalGrid.port_b,pipe1. port_a1) annotation (Line(points={{-60,-10},
              {-50,-10},{-50,-24},{-40,-24}}, color={0,127,255}));
      connect(externalGrid.port_a,pipe1. ports_b2[1]) annotation (Line(points={{-80,-10},
              {-90,-10},{-90,-36},{-40,-36}},      color={0,127,255}));
      connect(junction1.port_2,aggregatedConsumer. port_a) annotation (Line(points={{20,-24},
              {50,-24},{50,-10},{60,-10}},          color={0,127,255}));
      connect(aggregatedConsumer.port_b,junction2. port_1) annotation (Line(points={{80,-10},
              {90,-10},{90,-36},{40,-36}},          color={0,127,255}));
      connect(aggregatedConsumer.port_ht,soilTemperature. port) annotation (Line(
            points={{70,0},{70,30},{50,30},{50,40}},   color={191,0,0}));
      connect(pipe2.ports_b1[1], port_b) annotation (Line(points={{14,20},{14,
              80},{-50,80},{-50,100}},
                                   color={0,127,255}));
      connect(pipe2.port_a2, port_a) annotation (Line(points={{26,20},{26,80},{
              50,80},{50,100}},
                             color={0,127,255}));
      connect(port_a, port_a)
        annotation (Line(points={{50,100},{50,100}}, color={0,127,255}));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
            Text(
              extent={{-60,-110},{60,-170}},
              lineColor={0,0,0},
              pattern=LinePattern.None,
              lineThickness=1,
              fillColor={0,0,173},
              fillPattern=FillPattern.Solid,
              textString="%name"),         Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={0,0,0},
              fillColor={227,221,16},
              fillPattern=FillPattern.Solid,
              radius=20),
            Ellipse(
              extent={{-88,89},{88,-87}},
              lineColor={0,0,0},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),
            Text(
              extent={{-74,80},{86,-80}},
              lineColor={0,0,0},
              textString="P")}), Diagram(coordinateSystem(preserveAspectRatio=
                false)));
    end DHNetwork;

    model DHNetworkCtrl
      DisHeatLib.Boundary.SoilTemperature soilTemperature(
        inputType=DisHeatLib.Boundary.BaseClasses.InputTypeSoilTemp.Constant,
        T_const(displayUnit="degC") = 283.15,
        T_mean=283.15,
        T_amp(displayUnit="degC") = 10,
        t_min=0) annotation (Placement(transformation(extent={{40,60},{60,40}})));
      DisHeatLib.Pipes.DualPipe pipe1(
        show_T=true,
        redeclare package Medium = IBPSA.Media.Water,
        redeclare DisHeatLib.Pipes.Library.Isoplus.Isoplus_Std_DN100 pipeType,
        L=500,
        nPorts1=1,
        nPorts2=1)
        annotation (Placement(transformation(extent={{-40,-40},{-20,-20}})));
      DisHeatLib.Pipes.DualPipe pipe2(
        show_T=true,
        redeclare package Medium = IBPSA.Media.Water,
        redeclare DisHeatLib.Pipes.Library.Isoplus.Isoplus_Std_DN100 pipeType,
        L=10,
        nPorts2=1,
        nPorts1=1) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=90,
            origin={20,10})));
      IBPSA.Fluid.FixedResistances.Junction junction2(
        redeclare package Medium = IBPSA.Media.Water,
        m_flow_nominal={10,-20,10},
        dp_nominal={10,-10,10})
        annotation (Placement(transformation(extent={{40,-26},{20,-46}})));
      IBPSA.Fluid.FixedResistances.Junction junction1(
        redeclare package Medium = IBPSA.Media.Water,
        m_flow_nominal={20,-10,-10},
        dp_nominal={10,-10,-10})
        annotation (Placement(transformation(extent={{0,-14},{20,-34}})));
      Components.ExternalGridConstantPressure externalGrid(show_T=true)
        annotation (Placement(transformation(extent={{-80,-20},{-60,0}})));
      Components.AggregatedConsumers aggregatedConsumer(fileName=fileName,
        tableName=tableName,                            redeclare
          DisHeatLib.Pipes.Library.Isoplus.Isoplus_Std_DN50 pipeType)
        annotation (Placement(transformation(extent={{60,-20},{80,0}})));
      Modelica.Fluid.Interfaces.FluidPort_a port_a(redeclare package Medium =
            IBPSA.Media.Water)
        annotation (Placement(transformation(extent={{40,90},{60,110}})));
      Modelica.Fluid.Interfaces.FluidPort_b port_b(redeclare package Medium =
            IBPSA.Media.Water)
        annotation (Placement(transformation(extent={{-60,90},{-40,110}})));
      Modelica.Blocks.Sources.RealExpression m_flow_supply(y=externalGrid.m_flow)
        annotation (Placement(transformation(extent={{-80,-80},{-60,-60}})));
      Modelica.Blocks.Logical.Hysteresis hysteresis(uLow=m_flow_supply_low, uHigh=
            m_flow_supply_high)
        annotation (Placement(transformation(extent={{-40,-80},{-20,-60}})));
      Modelica.Blocks.Interfaces.BooleanOutput T_supply_secondary_ctrl annotation (
          Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=-90,
            origin={0,-110})));

      parameter Real m_flow_supply_low=2.5
        "if y=true and u<=uLow, switch to y=false";
      parameter Real m_flow_supply_high=3.0
        "if y=false and u>=uHigh, switch to y=true";

      parameter String fileName="modelica://TestbedExample/data/network_data.txt"
        "File where matrix is stored" annotation(Dialog(group="Demand profile"),loadSelector(filter="Text files(*.txt);;CSV files (*.csv)",caption="Open data file"));

      parameter String tableName="NetworkData" "Table name on file"  annotation(Dialog(group="Demand profile"));

    equation
      connect(soilTemperature.port,pipe1. port_ht) annotation (Line(points={{50,40},
              {50,30},{-30,30},{-30,-20}}, color={191,0,0}));
      connect(junction1.port_1,pipe1. ports_b1[1])
        annotation (Line(points={{0,-24},{-20,-24}}, color={0,127,255}));
      connect(junction1.port_3,pipe2. port_a1) annotation (Line(points={{10,-14},{10,
              -8},{14,-8},{14,0}},     color={0,127,255}));
      connect(pipe2.port_ht,soilTemperature. port) annotation (Line(points={{10,10},
              {-30,10},{-30,30},{50,30},{50,40}},  color={191,0,0}));
      connect(pipe1.port_a2,junction2. port_2)
        annotation (Line(points={{-20,-36},{20,-36}}, color={0,127,255}));
      connect(pipe2.ports_b2[1],junction2. port_3) annotation (Line(points={{26,0},{
              26,-8},{30,-8},{30,-26}},    color={0,127,255}));
      connect(externalGrid.port_b,pipe1. port_a1) annotation (Line(points={{-60,-10},
              {-50,-10},{-50,-24},{-40,-24}}, color={0,127,255}));
      connect(externalGrid.port_a,pipe1. ports_b2[1]) annotation (Line(points={{-80,-10},
              {-90,-10},{-90,-36},{-40,-36}},      color={0,127,255}));
      connect(junction1.port_2,aggregatedConsumer. port_a) annotation (Line(points={{20,-24},
              {50,-24},{50,-10},{60,-10}},          color={0,127,255}));
      connect(aggregatedConsumer.port_b,junction2. port_1) annotation (Line(points={{80,-10},
              {90,-10},{90,-36},{40,-36}},          color={0,127,255}));
      connect(aggregatedConsumer.port_ht,soilTemperature. port) annotation (Line(
            points={{70,0},{70,30},{50,30},{50,40}},   color={191,0,0}));
      connect(pipe2.ports_b1[1], port_b) annotation (Line(points={{14,20},{14,80},{-50,
              80},{-50,100}},      color={0,127,255}));
      connect(pipe2.port_a2, port_a) annotation (Line(points={{26,20},{26,80},{50,80},
              {50,100}},     color={0,127,255}));
      connect(port_a, port_a)
        annotation (Line(points={{50,100},{50,100}}, color={0,127,255}));
      connect(hysteresis.u, m_flow_supply.y)
        annotation (Line(points={{-42,-70},{-59,-70}}, color={0,0,127}));
      connect(hysteresis.y, T_supply_secondary_ctrl)
        annotation (Line(points={{-19,-70},{0,-70},{0,-110}}, color={255,0,255}));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
            Text(
              extent={{-60,-110},{60,-170}},
              lineColor={0,0,0},
              pattern=LinePattern.None,
              lineThickness=1,
              fillColor={0,0,173},
              fillPattern=FillPattern.Solid,
              textString="%name"),         Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={0,0,0},
              fillColor={227,221,16},
              fillPattern=FillPattern.Solid,
              radius=20),
            Ellipse(
              extent={{-88,89},{88,-87}},
              lineColor={0,0,0},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),
            Text(
              extent={{-74,80},{86,-80}},
              lineColor={0,0,0},
              textString="P")}), Diagram(coordinateSystem(preserveAspectRatio=
                false)));
    end DHNetworkCtrl;

    package Utility
      "Helper functions."
      function logMean "Logarithmic mean of two values"
        extends Modelica.Icons.Function;
        input Real a "First input value";
        input Real b "Second input value";
        output Real y "Result";

      protected
        Real delta = a - b;

      algorithm
          if abs(delta) > 1e-3 then
            y := delta / log(a/b);
          else
            y := a - delta/2*(1 + delta/6/a*(1 + delta/2/a));
          end if;

      end logMean;

      function clamp
        "Ensure that value x lies within the closed interval [a, b]"
        extends Modelica.Icons.Function;
        input Real a "Closed interval left limit";
        input Real x "Value limited by closed interval [a, b]";
        input Real b "Closed interval right limit";
        output Real y "Result";

      algorithm

        if x > a then
          if x < b then
            y := x;
          else
            y := b;
          end if;
        else
          y := a;
        end if;

      end clamp;

      function eta
        "Compute overall heat pump effciency"
        extends Modelica.Icons.Function;
        input Real T_evap_log "Evaporator log mean temperature";
        input Real T_cond_log "Condenser log mean temperature";
        input Real eta_sys "Relation between work provided by the pump and available thermodynamic work";
        output Real y "Result";

      algorithm
        if T_evap_log < T_cond_log then
          y := eta_sys/(1 - T_evap_log/T_cond_log);
        else
          y := Modelica.Constants.inf;
        end if;

      end eta;
      annotation (Icon(graphics={
            Rectangle(
              lineColor={200,200,200},
              fillColor={248,248,248},
              fillPattern=FillPattern.HorizontalCylinder,
              extent={{-100,-100},{100,100}},
              radius=25),                  Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={0,0,0},
              radius=25),
        Polygon(
          origin={-2.6165,-0.1418},
          rotation=45.0,
          fillColor={64,64,64},
          pattern=LinePattern.None,
          fillPattern=FillPattern.Solid,
          points={{-15.0,93.333},{-15.0,68.333},{0.0,58.333},{15.0,68.333},{15.0,93.333},{20.0,93.333},{25.0,83.333},{25.0,58.333},{10.0,43.333},{10.0,-41.667},{25.0,-56.667},{25.0,-76.667},{10.0,-91.667},{0.0,-91.667},{0.0,-81.667},{5.0,-81.667},{15.0,-71.667},{15.0,-61.667},{5.0,-51.667},{-5.0,-51.667},{-15.0,-61.667},{-15.0,-71.667},{-5.0,-81.667},{0.0,-81.667},{0.0,-91.667},{-10.0,-91.667},{-25.0,-76.667},{-25.0,-56.667},{-10.0,-41.667},{-10.0,43.333},{-25.0,58.333},{-25.0,83.333},{-20.0,93.333}}),
        Polygon(
          origin={6.1018,9.218},
          rotation=-45.0,
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid,
          points={{-15.0,87.273},{15.0,87.273},{20.0,82.273},{20.0,27.273},{10.0,17.273},{10.0,7.273},{20.0,2.273},{20.0,-2.727},{5.0,-2.727},{5.0,-77.727},{10.0,-87.727},{5.0,-112.727},{-5.0,-112.727},{-10.0,-87.727},{-5.0,-77.727},{-5.0,-2.727},{-20.0,-2.727},{-20.0,2.273},{-10.0,7.273},{-10.0,17.273},{-20.0,27.273},{-20.0,82.273}})}));
    end Utility;
    annotation (Icon(graphics={
          Rectangle(
            lineColor={200,200,200},
            fillColor={248,248,248},
            fillPattern=FillPattern.HorizontalCylinder,
            extent={{-100,-100},{100,100}},
            radius=25),                  Rectangle(
            extent={{-100,100},{100,-100}},
            lineColor={0,0,0},
            radius=25),
          Bitmap(extent={{-64,-76},{76,64}}, fileName=
                "modelica://TestbedExample/images/gears.svg")}));
  end Components;

  package TestCases

    package TestCase1
      model TestCase1_ReferenceSystem

        parameter Modelica.SIunits.MassFlowRate m1_flow_nominal = 5 "Nominal mass flow rate primary side"
          annotation(Dialog(group="Nominal condition"));
        parameter Modelica.SIunits.MassFlowRate m2_flow_nominal = 5 "Nominal mass flow rate secondary side"
          annotation(Dialog(group="Nominal condition"));
        parameter Modelica.SIunits.PressureDifference dp1_nominal(displayUnit=
              "bar") = 200000 "Nominal pressure difference primary side"
              annotation(Dialog(group="Nominal condition"));
        parameter Modelica.SIunits.PressureDifference dp2_nominal(displayUnit=
              "bar") = 200000 "Nominal pressure difference secondary side"
              annotation(Dialog(group="Nominal condition"));
        parameter Modelica.SIunits.PressureDifference dpValve_nominal=200000
          "Nominal pressure drop of fully open valve" annotation(Dialog(group="Nominal condition"));

        parameter Modelica.SIunits.Power Q_flow_nominal(displayUnit="kW") = 50000
          "Nominal heat flow rate" annotation(Dialog(group="Nominal condition"));
        parameter Modelica.SIunits.Temperature T_supply_primary_nominal=348.15
          "Nominal primary supply temperature" annotation(Dialog(group="Nominal condition"));
        parameter Modelica.SIunits.Temperature T_return_primary_nominal=308.15
          "Nominal primary return temperature" annotation(Dialog(group="Nominal condition"));
        parameter Modelica.SIunits.Temperature T_supply_secondary_nominal=341.15 "Nominal secondary supply temperature" annotation(Dialog(group="Nominal condition"));
        parameter Modelica.SIunits.Temperature T_return_secondary_nominal=301.15 "Nominal secondary return temperature" annotation(Dialog(group="Nominal condition"));

        DisHeatLib.Boundary.SoilTemperature soilTemperature(
          inputType=DisHeatLib.Boundary.BaseClasses.InputTypeSoilTemp.Constant,
          T_const(displayUnit="degC") = 283.15,
          T_mean=283.15,
          T_amp(displayUnit="degC") = 10,
          t_min=0) annotation (Placement(transformation(extent={{10,90},{30,70}})));
        DisHeatLib.Pipes.DualPipe pipe1(
          show_T=true,
          redeclare package Medium = IBPSA.Media.Water,
          redeclare DisHeatLib.Pipes.Library.Isoplus.Isoplus_Std_DN100 pipeType,
          L=500,
          nPorts1=1,
          nPorts2=1)
          annotation (Placement(transformation(extent={{-40,-80},{-20,-60}})));
        DisHeatLib.Pipes.DualPipe pipe2(
          show_T=true,
          redeclare package Medium = IBPSA.Media.Water,
          redeclare DisHeatLib.Pipes.Library.Isoplus.Isoplus_Std_DN100 pipeType,
          L=10,
          nPorts2=1,
          nPorts1=1) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=90,
              origin={20,-30})));
        IBPSA.Fluid.FixedResistances.Junction junction2(
          redeclare package Medium = IBPSA.Media.Water,
          m_flow_nominal={10,-20,10},
          dp_nominal={10,-10,10})
          annotation (Placement(transformation(extent={{40,-66},{20,-86}})));
        IBPSA.Fluid.FixedResistances.Junction junction1(
          redeclare package Medium = IBPSA.Media.Water,
          m_flow_nominal={20,-10,-10},
          dp_nominal={10,-10,-10})
          annotation (Placement(transformation(extent={{0,-54},{20,-74}})));
        Components.ExternalGridConstantPressure externalGrid(show_T=
              true)
          annotation (Placement(transformation(extent={{-80,-60},{-60,-40}})));
        Components.AggregatedConsumers aggregatedConsumer(redeclare
            DisHeatLib.Pipes.Library.Isoplus.Isoplus_Std_DN50 pipeType)
          annotation (Placement(transformation(extent={{60,-60},{80,-40}})));
        Components.SingleConsumer singleConsumer(
          m_flow_nominal=m2_flow_nominal,
          dp_nominal=dp2_nominal,
          T_supply_nominal=T_supply_secondary_nominal,
          T_return_nominal=T_return_secondary_nominal,
          Q_flow_nominal=Q_flow_nominal,
          scaling=2,
          fan(show_T=true, addPowerToMedium=false))
          annotation (Placement(transformation(extent={{10,30},{30,50}})));
        SubstationTeststand.Components.Substation.Substation substation(
          m1_flow_nominal=m1_flow_nominal,
          m2_flow_nominal=m2_flow_nominal,
          dp1_nominal=dp1_nominal,
          dp2_nominal=dp2_nominal,
          dpValve_nominal=dpValve_nominal,
          hex(show_T=true))
          annotation (Placement(transformation(extent={{10,20},{30,0}})));

      equation

        substation.T_supply_secondary_set = T_supply_secondary_nominal;
        connect(soilTemperature.port, pipe1.port_ht) annotation (Line(points={{20,70},
                {20,60},{-30,60},{-30,-60}}, color={191,0,0}));
        connect(junction1.port_1, pipe1.ports_b1[1])
          annotation (Line(points={{0,-64},{-20,-64}}, color={0,127,255}));
        connect(junction1.port_3, pipe2.port_a1) annotation (Line(points={{10,-54},{10,
                -48},{14,-48},{14,-40}}, color={0,127,255}));
        connect(pipe2.port_ht, soilTemperature.port) annotation (Line(points={{10,-30},
                {-30,-30},{-30,60},{20,60},{20,70}}, color={191,0,0}));
        connect(pipe1.port_a2, junction2.port_2)
          annotation (Line(points={{-20,-76},{20,-76}}, color={0,127,255}));
        connect(pipe2.ports_b2[1], junction2.port_3) annotation (Line(points={{26,-40},
                {26,-48},{30,-48},{30,-66}}, color={0,127,255}));
        connect(externalGrid.port_b, pipe1.port_a1) annotation (Line(points={{-60,-50},
                {-50,-50},{-50,-64},{-40,-64}}, color={0,127,255}));
        connect(externalGrid.port_a, pipe1.ports_b2[1]) annotation (Line(points={{-80,-50},
                {-90,-50},{-90,-76},{-40,-76}},      color={0,127,255}));
        connect(junction1.port_2, aggregatedConsumer.port_a) annotation (Line(points={{20,-64},
                {50,-64},{50,-50},{60,-50}},          color={0,127,255}));
        connect(aggregatedConsumer.port_b, junction2.port_1) annotation (Line(points={{80,-50},
                {90,-50},{90,-76},{40,-76}},          color={0,127,255}));
        connect(aggregatedConsumer.port_ht, soilTemperature.port) annotation (Line(
              points={{70,-40},{70,60},{20,60},{20,70}}, color={191,0,0}));
        connect(substation.port_a2, singleConsumer.port_b) annotation (Line(points={{30,
                16},{40,16},{40,40},{30,40}}, color={0,127,255}));
        connect(substation.port_b2, singleConsumer.port_a) annotation (Line(points={{10,
                16},{0,16},{0,40},{10,40}}, color={0,127,255}));
        connect(substation.port_b1, pipe2.port_a2) annotation (Line(points={{30,4},{40,
                4},{40,-10},{26,-10},{26,-20}}, color={0,127,255}));
        connect(substation.port_a1, pipe2.ports_b1[1]) annotation (Line(points={{10,4},
                {0,4},{0,-10},{14,-10},{14,-20}}, color={0,127,255}));
        annotation (
          Icon(
            coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
                  100}}), graphics={
              Ellipse(lineColor = {75,138,73},
                      fillColor={255,255,255},
                      fillPattern = FillPattern.Solid,
                      extent={{-100,-100},{100,100}}),
              Polygon(lineColor = {0,0,255},
                      fillColor = {75,138,73},
                      pattern = LinePattern.None,
                      fillPattern = FillPattern.Solid,
                      points={{-36,60},{64,0},{-36,-60},{-36,60}})}),
          Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                  {100,100}})),
          experiment(
            StopTime=172800,
            Interval=10,
            __Dymola_Algorithm="Dassl"),
          Documentation(info="<html>
    <p>This model implements the thermal system of the JRA 1.1 multi-energy benchmark. The heating network is operated by controlling the pressure, keeping the pressure drop at the consumers above a certain threshold.</p>
    </html>"));
      end TestCase1_ReferenceSystem;

      model TestCase1_ReferenceSystem_subsystems

        parameter Modelica.SIunits.MassFlowRate m1_flow_nominal = 5 "Nominal mass flow rate primary side"
          annotation(Dialog(group="Nominal condition"));
        parameter Modelica.SIunits.MassFlowRate m2_flow_nominal = 5 "Nominal mass flow rate secondary side"
          annotation(Dialog(group="Nominal condition"));
        parameter Modelica.SIunits.PressureDifference dp1_nominal(displayUnit=
              "bar") = 200000 "Nominal pressure difference primary side"
              annotation(Dialog(group="Nominal condition"));
        parameter Modelica.SIunits.PressureDifference dp2_nominal(displayUnit=
              "bar") = 200000 "Nominal pressure difference secondary side"
              annotation(Dialog(group="Nominal condition"));
        parameter Modelica.SIunits.PressureDifference dpValve_nominal=200000
          "Nominal pressure drop of fully open valve" annotation(Dialog(group="Nominal condition"));

        parameter Modelica.SIunits.Power Q_flow_nominal(displayUnit="kW") = 50000
          "Nominal heat flow rate" annotation(Dialog(group="Nominal condition"));
        parameter Modelica.SIunits.Temperature T_supply_primary_nominal=348.15
          "Nominal primary supply temperature" annotation(Dialog(group="Nominal condition"));
        parameter Modelica.SIunits.Temperature T_return_primary_nominal=308.15
          "Nominal primary return temperature" annotation(Dialog(group="Nominal condition"));
        parameter Modelica.SIunits.Temperature T_supply_secondary_nominal=341.15 "Nominal secondary supply temperature" annotation(Dialog(group="Nominal condition"));
        parameter Modelica.SIunits.Temperature T_return_secondary_nominal=301.15 "Nominal secondary return temperature" annotation(Dialog(group="Nominal condition"));

        Components.DHNetwork dHNetwork
          annotation (Placement(transformation(extent={{-10,-60},{10,-40}})));
        SubstationTeststand.Components.Substation.Substation substation(
          m1_flow_nominal=m1_flow_nominal,
          m2_flow_nominal=m2_flow_nominal,
          dp1_nominal=dp1_nominal,
          dp2_nominal=dp2_nominal,
          dpValve_nominal=dpValve_nominal,
          hex(show_T=true))
          annotation (Placement(transformation(extent={{-10,20},{10,0}})));
        Components.SingleConsumer singleConsumer(
          m_flow_nominal=m2_flow_nominal,
          dp_nominal=dp2_nominal,
          T_supply_nominal=T_supply_secondary_nominal,
          T_return_nominal=T_return_secondary_nominal,
          Q_flow_nominal=Q_flow_nominal,
          scaling=2,
          fan(show_T=true, addPowerToMedium=false))
          annotation (Placement(transformation(extent={{-10,40},{10,60}})));

      equation

        substation.T_supply_secondary_set = T_supply_secondary_nominal;

        connect(substation.port_b2, singleConsumer.port_a) annotation (Line(points={{-10,
                16},{-30,16},{-30,50},{-10,50}}, color={0,127,255}));
        connect(singleConsumer.port_b, substation.port_a2) annotation (Line(points={{10,
                50},{30,50},{30,16},{10,16}}, color={0,127,255}));
        connect(substation.port_b1, dHNetwork.port_a) annotation (Line(points={{10,4},
                {30,4},{30,-30},{5,-30},{5,-40}}, color={0,127,255}));
        connect(substation.port_a1, dHNetwork.port_b) annotation (Line(points={{-10,4},
                {-30,4},{-30,-30},{-5,-30},{-5,-40}}, color={0,127,255}));
        annotation (
          Icon(
            coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
                  100}}), graphics={
              Ellipse(lineColor = {75,138,73},
                      fillColor={255,255,255},
                      fillPattern = FillPattern.Solid,
                      extent={{-100,-100},{100,100}}),
              Polygon(lineColor = {0,0,255},
                      fillColor = {75,138,73},
                      pattern = LinePattern.None,
                      fillPattern = FillPattern.Solid,
                      points={{-36,60},{64,0},{-36,-60},{-36,60}})}),
          Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                  {100,100}})),
          experiment(
            StopTime=172800,
            __Dymola_NumberOfIntervals=2880,
            __Dymola_Algorithm="Dassl"),
          Documentation(info="<html>
    <p>This model implements the thermal system of the JRA 1.1 multi-energy benchmark. The heating network is operated by controlling the pressure, keeping the pressure drop at the consumers above a certain threshold.</p>
    </html>"));
      end TestCase1_ReferenceSystem_subsystems;

      model TestCase1_Testbed

        parameter Modelica.SIunits.MassFlowRate m1_flow_nominal = 5 "Nominal mass flow rate primary side"
          annotation(Dialog(group="Nominal condition"));
        parameter Modelica.SIunits.MassFlowRate m2_flow_nominal = 5 "Nominal mass flow rate secondary side"
          annotation(Dialog(group="Nominal condition"));
        parameter Modelica.SIunits.PressureDifference dp1_nominal(displayUnit=
              "bar") = 200000 "Nominal pressure difference primary side"
              annotation(Dialog(group="Nominal condition"));
        parameter Modelica.SIunits.PressureDifference dp2_nominal(displayUnit=
              "bar") = 200000 "Nominal pressure difference secondary side"
              annotation(Dialog(group="Nominal condition"));
        parameter Modelica.SIunits.PressureDifference dpValve_nominal=200000
          "Nominal pressure drop of fully open valve" annotation(Dialog(group="Nominal condition"));

        parameter Modelica.SIunits.PressureDifference delta_p_secondary_nominal(displayUnit="bar") = 200000 "Nominal pressure difference for secondary side";
        parameter Modelica.SIunits.PressureDifference p_ref_secondary_nominal(displayUnit="bar") = IBPSA.Media.Water.p_default "Reference pressure level for secondary side";

        parameter Modelica.SIunits.Power Q_flow_nominal(displayUnit="kW") = 50000
          "Nominal heat flow rate" annotation(Dialog(group="Nominal condition"));
        parameter Modelica.SIunits.Temperature T_supply_primary_nominal=348.15
          "Nominal primary supply temperature" annotation(Dialog(group="Nominal condition"));
        parameter Modelica.SIunits.Temperature T_return_primary_nominal=308.15
          "Nominal primary return temperature" annotation(Dialog(group="Nominal condition"));
        parameter Modelica.SIunits.Temperature T_supply_secondary_nominal=341.15 "Nominal secondary supply temperature" annotation(Dialog(group="Nominal condition"));
        parameter Modelica.SIunits.Temperature T_return_secondary_nominal=301.15 "Nominal secondary return temperature" annotation(Dialog(group="Nominal condition"));

        Components.DHNetwork dHNetwork
          annotation (Placement(transformation(extent={{-10,-90},{10,-70}})));
        DisHeatLib.Utilities.ThermoHydraulicEquivalent.Demand_mT demand_mT(
            m_flow_nominal=4, dp_nominal=200000)
          annotation (Placement(transformation(extent={{-6,-46},{6,-34}})));
        DisHeatLib.Utilities.ThermoHydraulicEquivalent.Supply_pT supply_pT
          annotation (Placement(transformation(extent={{-6,24},{6,36}})));
        SubstationTeststand.TestStandSimplifiedPIDPressureCtrl
                                           testStandSimplifiedPIDPressureCtrl(
          m1_flow_nominal=m1_flow_nominal,
          m2_flow_nominal=m2_flow_nominal,
          dp1_nominal=dp1_nominal,
          dp2_nominal=delta_p_secondary_nominal,
          dpValve_nominal=dpValve_nominal,
          substation(hex(show_T=true)))
          annotation (Placement(transformation(extent={{-10,0},{10,-20}})));
        Modelica.Blocks.Sources.Constant T_supply_secondary(k=
              T_supply_secondary_nominal)
          annotation (Placement(transformation(extent={{-70,-20},{-50,0}})));
        Modelica.Blocks.Sources.Constant delta_p_secondary(k=
              delta_p_secondary_nominal)
          annotation (Placement(transformation(extent={{70,-20},{50,0}})));
        Components.SingleConsumer singleConsumer(
          use_p_ref=false,
          m_flow_nominal=m2_flow_nominal,
          dp_nominal=delta_p_secondary_nominal,
          T_supply_nominal=T_supply_secondary_nominal,
          T_return_nominal=T_return_secondary_nominal,
          Q_flow_nominal=Q_flow_nominal,
          scaling=2,
          fan(show_T=true, addPowerToMedium=false))
          annotation (Placement(transformation(extent={{10,70},{-10,90}})));
      equation

        supply_pT.delta_p_set = 0;
        supply_pT.p_ref_set = p_ref_secondary_nominal;

        connect(demand_mT.port_b, dHNetwork.port_a) annotation (Line(points={{6,
                -40},{20,-40},{20,-60},{5,-60},{5,-70}}, color={0,127,255}));
        connect(demand_mT.port_a, dHNetwork.port_b) annotation (Line(points={{-6,
                -40},{-20,-40},{-20,-60},{-5,-60},{-5,-70}}, color={0,127,255}));
        connect(demand_mT.T_supply, testStandSimplifiedPIDPressureCtrl.T_supply_primary_set)
          annotation (Line(points={{-3.6,-46.6},{-3.6,-50},{-30,-50},{-30,-16},
                {-12,-16}}, color={0,0,127}));
        connect(supply_pT.T_return, testStandSimplifiedPIDPressureCtrl.T_return_secondary_set)
          annotation (Line(points={{-2.4,23.4},{-2.4,10},{-30,10},{-30,-4},{-12,
                -4}}, color={0,0,127}));
        connect(testStandSimplifiedPIDPressureCtrl.m_flow_return_secondary,
          supply_pT.m_flow) annotation (Line(points={{12,-10},{30,-10},{30,10},
                {2.4,10},{2.4,23.4}}, color={0,0,127}));
        connect(testStandSimplifiedPIDPressureCtrl.m_flow_return_primary,
          demand_mT.m_flow_set) annotation (Line(points={{-6,-21},{-6,-28},{-2.4,
                -28},{-2.4,-32.8}}, color={0,0,127}));
        connect(testStandSimplifiedPIDPressureCtrl.T_return_primary, demand_mT.T_return_set)
          annotation (Line(points={{6,-21},{6,-28},{2.4,-28},{2.4,-32.8}},
              color={0,0,127}));
        connect(testStandSimplifiedPIDPressureCtrl.T_supply_secondary,
          supply_pT.T_supply_set) annotation (Line(points={{6,1},{6,16},{12,16},
                {12,46},{3.6,46},{3.6,37.2}}, color={0,0,127}));
        connect(testStandSimplifiedPIDPressureCtrl.delta_p_secondary_set,
          delta_p_secondary.y) annotation (Line(points={{12,-4},{40,-4},{40,-10},
                {49,-10}}, color={0,0,127}));
        connect(T_supply_secondary.y, testStandSimplifiedPIDPressureCtrl.T_supply_secondary_set)
          annotation (Line(points={{-49,-10},{-12,-10}}, color={0,0,127}));
        connect(testStandSimplifiedPIDPressureCtrl.delta_p_primary_set,
          demand_mT.delta_p) annotation (Line(points={{12,-16},{30,-16},{30,-50},
                {3.6,-50},{3.6,-46.6}}, color={0,0,127}));
        connect(singleConsumer.port_a, supply_pT.port_b) annotation (Line(points={{10,
                80},{30,80},{30,30},{6,30}}, color={0,127,255}));
        connect(singleConsumer.port_b, supply_pT.port_a) annotation (Line(points={{-10,
                80},{-32,80},{-32,30},{-6,30}}, color={0,127,255}));
        annotation (
          Icon(
            coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
                  100}}), graphics={
              Ellipse(lineColor = {75,138,73},
                      fillColor={255,255,255},
                      fillPattern = FillPattern.Solid,
                      extent={{-100,-100},{100,100}}),
              Polygon(lineColor = {0,0,255},
                      fillColor = {75,138,73},
                      pattern = LinePattern.None,
                      fillPattern = FillPattern.Solid,
                      points={{-36,60},{64,0},{-36,-60},{-36,60}})}),
          Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                  {100,100}})),
          experiment(
            StopTime=172800,
            __Dymola_NumberOfIntervals=2880,
            __Dymola_Algorithm="Dassl"),
          Documentation(info="<html>
    <p>This model implements the thermal system of the JRA 1.1 multi-energy benchmark. The heating network is operated by controlling the pressure, keeping the pressure drop at the consumers above a certain threshold.</p>
    </html>"));
      end TestCase1_Testbed;

      model TestCase1_Testbed_subsystems
        FMI.Teststand teststand
          annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
        FMI.SingleConsumer consumer
          annotation (Placement(transformation(extent={{-10,50},{10,70}})));
        FMI.DHNetwork network
          annotation (Placement(transformation(extent={{-10,-70},{10,-50}})));
      equation
        connect(teststand.T_supply_secondary, consumer.T_supply_set) annotation (Line(
              points={{0,11},{0,20},{40,20},{40,90},{0,90},{0,72}}, color={0,0,127}));
        connect(consumer.T_return, teststand.T_return_secondary_set) annotation (Line(
              points={{-6,49},{-6,40},{-40,40},{-40,6},{-12,6}}, color={0,0,127}));
        connect(consumer.delta_p_secondary, teststand.delta_p_secondary_set)
          annotation (Line(points={{0,49},{0,34},{20,34},{20,5.8},{12,5.8}}, color={0,
                0,127}));
        connect(consumer.m_flow, teststand.m_flow_return_secondary) annotation (Line(
              points={{6,49},{6,40},{26,40},{26,0},{12,0}}, color={0,0,127}));
        connect(teststand.T_return_primary, network.T_return_set)
          annotation (Line(points={{4,-11},{4,-48}}, color={0,0,127}));
        connect(teststand.m_flow_return_primary, network.m_flow) annotation (
            Line(points={{-4,-11},{-4,-30},{-4,-48},{-3.8,-48}}, color={0,0,127}));
        connect(network.delta_p, teststand.delta_p_primary_set) annotation (
            Line(points={{4,-71},{4,-88},{40,-88},{40,-6},{12,-6}}, color={0,0,127}));
        connect(network.T_supply, teststand.T_supply_primary_set) annotation (
           Line(points={{-3.8,-71},{-3.8,-88},{-40,-88},{-40,-6},{-12,-6}}, color={0,0,
                127}));
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
          experiment(
            StopTime=172800,
            __Dymola_NumberOfIntervals=2880,
            __Dymola_Algorithm="Dassl"));
      end TestCase1_Testbed_subsystems;

    end TestCase1;

    package TestCase2
      model TestCase2_ReferenceSystem

        parameter Modelica.SIunits.MassFlowRate m1_flow_nominal = 5 "Nominal mass flow rate primary side"
          annotation(Dialog(group="Nominal condition"));
        parameter Modelica.SIunits.MassFlowRate m2_flow_nominal = 5 "Nominal mass flow rate secondary side"
          annotation(Dialog(group="Nominal condition"));
        parameter Modelica.SIunits.PressureDifference dp1_nominal(displayUnit=
              "bar") = 200000 "Nominal pressure difference primary side"
              annotation(Dialog(group="Nominal condition"));
        parameter Modelica.SIunits.PressureDifference dp2_nominal(displayUnit=
              "bar") = 200000 "Nominal pressure difference secondary side"
              annotation(Dialog(group="Nominal condition"));
        parameter Modelica.SIunits.PressureDifference dpValve_nominal=200000
          "Nominal pressure drop of fully open valve" annotation(Dialog(group="Nominal condition"));

        parameter Modelica.SIunits.Power Q_flow_nominal(displayUnit="kW") = 50000
          "Nominal heat flow rate" annotation(Dialog(group="Nominal condition"));
        parameter Modelica.SIunits.Temperature T_supply_primary_nominal=348.15
          "Nominal primary supply temperature" annotation(Dialog(group="Nominal condition"));
        parameter Modelica.SIunits.Temperature T_return_primary_nominal=308.15
          "Nominal primary return temperature" annotation(Dialog(group="Nominal condition"));
        parameter Modelica.SIunits.Temperature T_supply_secondary_nominal=341.15 "Nominal secondary supply temperature" annotation(Dialog(group="Nominal condition"));
        parameter Modelica.SIunits.Temperature T_return_secondary_nominal=301.15 "Nominal secondary return temperature" annotation(Dialog(group="Nominal condition"));

        DisHeatLib.Boundary.SoilTemperature soilTemperature(
          inputType=DisHeatLib.Boundary.BaseClasses.InputTypeSoilTemp.Constant,
          T_const(displayUnit="degC") = 283.15,
          T_mean=283.15,
          T_amp(displayUnit="degC") = 10,
          t_min=0) annotation (Placement(transformation(extent={{10,90},{30,70}})));
        DisHeatLib.Pipes.DualPipe pipe1(
          show_T=true,
          redeclare package Medium = IBPSA.Media.Water,
          redeclare DisHeatLib.Pipes.Library.Isoplus.Isoplus_Std_DN100 pipeType,
          L=500,
          nPorts1=1,
          nPorts2=1)
          annotation (Placement(transformation(extent={{-40,-80},{-20,-60}})));
        DisHeatLib.Pipes.DualPipe pipe2(
          show_T=true,
          redeclare package Medium = IBPSA.Media.Water,
          redeclare DisHeatLib.Pipes.Library.Isoplus.Isoplus_Std_DN100 pipeType,
          L=10,
          nPorts2=1,
          nPorts1=1) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=90,
              origin={20,-30})));
        IBPSA.Fluid.FixedResistances.Junction junction2(
          redeclare package Medium = IBPSA.Media.Water,
          m_flow_nominal={10,-20,10},
          dp_nominal={10,-10,10})
          annotation (Placement(transformation(extent={{40,-66},{20,-86}})));
        IBPSA.Fluid.FixedResistances.Junction junction1(
          redeclare package Medium = IBPSA.Media.Water,
          m_flow_nominal={20,-10,-10},
          dp_nominal={10,-10,-10})
          annotation (Placement(transformation(extent={{0,-54},{20,-74}})));
        Components.ExternalGridConstantPressure externalGrid(show_T=
              true)
          annotation (Placement(transformation(extent={{-80,-60},{-60,-40}})));
        Components.AggregatedConsumers aggregatedConsumer(redeclare
            DisHeatLib.Pipes.Library.Isoplus.Isoplus_Std_DN50 pipeType)
          annotation (Placement(transformation(extent={{60,-60},{80,-40}})));
        Components.SingleConsumerWithBooster
                                  singleConsumerWithBooster(
          m_flow_nominal=m2_flow_nominal,
          dp_nominal=dp2_nominal,
          T_supply_nominal=T_supply_secondary_nominal,
          T_return_nominal=T_return_secondary_nominal,
          Q_flow_nominal=Q_flow_nominal,
          scaling=2,
          fan(show_T=true, addPowerToMedium=false),
          heatpump(show_T=true, evaporator(dp_nominal(displayUnit="Pa"))))
          annotation (Placement(transformation(extent={{10,30},{30,50}})));
        SubstationTeststand.Components.Substation.Substation substation(
          m1_flow_nominal=m1_flow_nominal,
          m2_flow_nominal=m2_flow_nominal,
          dp1_nominal=dp1_nominal,
          dp2_nominal=dp2_nominal,
          dpValve_nominal=dpValve_nominal,
          hex(show_T=true))
          annotation (Placement(transformation(extent={{10,20},{30,0}})));

        Modelica.Blocks.Logical.Switch switch
          annotation (Placement(transformation(extent={{-14,6},{-6,14}})));
        Modelica.Blocks.Sources.RealExpression m_flow_supply(y=externalGrid.m_flow)
          annotation (Placement(transformation(extent={{-90,0},{-70,20}})));
        Modelica.Blocks.Sources.Constant T_supply_secondary_high(k=
              T_supply_secondary_nominal) annotation (Placement(transformation(
              extent={{-4,4},{4,-4}},
              rotation=90,
              origin={-20,-4})));
        Modelica.Blocks.Sources.Constant T_supply_secondary_low(k=
              T_supply_secondary_nominal - 5.0) annotation (Placement(
              transformation(
              extent={{-4,-4},{4,4}},
              rotation=-90,
              origin={-20,24})));
        Modelica.Blocks.Logical.Hysteresis hysteresis(uLow=2.5, uHigh=3.0)
          annotation (Placement(transformation(extent={{-60,0},{-40,20}})));
      equation

        connect(soilTemperature.port, pipe1.port_ht) annotation (Line(points={{20,70},
                {20,60},{-30,60},{-30,-60}}, color={191,0,0}));
        connect(junction1.port_1, pipe1.ports_b1[1])
          annotation (Line(points={{0,-64},{-20,-64}}, color={0,127,255}));
        connect(junction1.port_3, pipe2.port_a1) annotation (Line(points={{10,-54},{10,
                -48},{14,-48},{14,-40}}, color={0,127,255}));
        connect(pipe2.port_ht, soilTemperature.port) annotation (Line(points={{10,-30},
                {-30,-30},{-30,60},{20,60},{20,70}}, color={191,0,0}));
        connect(pipe1.port_a2, junction2.port_2)
          annotation (Line(points={{-20,-76},{20,-76}}, color={0,127,255}));
        connect(pipe2.ports_b2[1], junction2.port_3) annotation (Line(points={{26,-40},
                {26,-48},{30,-48},{30,-66}}, color={0,127,255}));
        connect(externalGrid.port_b, pipe1.port_a1) annotation (Line(points={{-60,-50},
                {-50,-50},{-50,-64},{-40,-64}}, color={0,127,255}));
        connect(externalGrid.port_a, pipe1.ports_b2[1]) annotation (Line(points={{-80,-50},
                {-90,-50},{-90,-76},{-40,-76}},      color={0,127,255}));
        connect(junction1.port_2, aggregatedConsumer.port_a) annotation (Line(points={{20,-64},
                {50,-64},{50,-50},{60,-50}},          color={0,127,255}));
        connect(aggregatedConsumer.port_b, junction2.port_1) annotation (Line(points={{80,-50},
                {90,-50},{90,-76},{40,-76}},          color={0,127,255}));
        connect(aggregatedConsumer.port_ht, soilTemperature.port) annotation (Line(
              points={{70,-40},{70,60},{20,60},{20,70}}, color={191,0,0}));
        connect(substation.port_a2, singleConsumerWithBooster.port_b) annotation (
            Line(points={{30,16},{40,16},{40,40},{30,40}}, color={0,127,255}));
        connect(substation.port_b2, singleConsumerWithBooster.port_a) annotation (
            Line(points={{10,16},{0,16},{0,40},{10,40}}, color={0,127,255}));
        connect(substation.port_b1, pipe2.port_a2) annotation (Line(points={{30,4},{40,
                4},{40,-10},{26,-10},{26,-20}}, color={0,127,255}));
        connect(substation.port_a1, pipe2.ports_b1[1]) annotation (Line(points={{10,4},
                {0,4},{0,-10},{14,-10},{14,-20}}, color={0,127,255}));
        connect(switch.y, substation.T_supply_secondary_set)
          annotation (Line(points={{-5.6,10},{8,10}}, color={0,0,127}));
        connect(hysteresis.u, m_flow_supply.y)
          annotation (Line(points={{-62,10},{-69,10}}, color={0,0,127}));
        connect(hysteresis.y, switch.u2)
          annotation (Line(points={{-39,10},{-14.8,10}}, color={255,0,255}));
        connect(T_supply_secondary_low.y, switch.u1) annotation (Line(points={{
                -20,19.6},{-20,13.2},{-14.8,13.2}}, color={0,0,127}));
        connect(T_supply_secondary_high.y, switch.u3) annotation (Line(points={
                {-20,0.4},{-20,6.8},{-14.8,6.8}}, color={0,0,127}));
        annotation (
          Icon(
            coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
                  100}}), graphics={
              Ellipse(lineColor = {75,138,73},
                      fillColor={255,255,255},
                      fillPattern = FillPattern.Solid,
                      extent={{-100,-100},{100,100}}),
              Polygon(lineColor = {0,0,255},
                      fillColor = {75,138,73},
                      pattern = LinePattern.None,
                      fillPattern = FillPattern.Solid,
                      points={{-36,60},{64,0},{-36,-60},{-36,60}})}),
          Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                  {100,100}})),
          experiment(
            StopTime=172800,
            Interval=10,
            __Dymola_Algorithm="Dassl"),
          Documentation(info="<html>
    <p>This model implements the thermal system of the JRA 1.1 multi-energy benchmark. The heating network is operated by controlling the pressure, keeping the pressure drop at the consumers above a certain threshold.</p>
    </html>"));
      end TestCase2_ReferenceSystem;

      model TestCase2_ReferenceSystem_subsystems

        parameter Modelica.SIunits.MassFlowRate m1_flow_nominal = 5 "Nominal mass flow rate primary side"
          annotation(Dialog(group="Nominal condition"));
        parameter Modelica.SIunits.MassFlowRate m2_flow_nominal = 5 "Nominal mass flow rate secondary side"
          annotation(Dialog(group="Nominal condition"));
        parameter Modelica.SIunits.PressureDifference dp1_nominal(displayUnit=
              "bar") = 200000 "Nominal pressure difference primary side"
              annotation(Dialog(group="Nominal condition"));
        parameter Modelica.SIunits.PressureDifference dp2_nominal(displayUnit=
              "bar") = 200000 "Nominal pressure difference secondary side"
              annotation(Dialog(group="Nominal condition"));
        parameter Modelica.SIunits.PressureDifference dpValve_nominal=200000
          "Nominal pressure drop of fully open valve" annotation(Dialog(group="Nominal condition"));

        parameter Modelica.SIunits.Power Q_flow_nominal(displayUnit="kW") = 50000
          "Nominal heat flow rate" annotation(Dialog(group="Nominal condition"));
        parameter Modelica.SIunits.Temperature T_supply_primary_nominal=348.15
          "Nominal primary supply temperature" annotation(Dialog(group="Nominal condition"));
        parameter Modelica.SIunits.Temperature T_return_primary_nominal=308.15
          "Nominal primary return temperature" annotation(Dialog(group="Nominal condition"));
        parameter Modelica.SIunits.Temperature T_supply_secondary_nominal=341.15 "Nominal secondary supply temperature" annotation(Dialog(group="Nominal condition"));
        parameter Modelica.SIunits.Temperature T_return_secondary_nominal=301.15 "Nominal secondary return temperature" annotation(Dialog(group="Nominal condition"));

        Components.DHNetworkCtrl dHNetworkCtrl
          annotation (Placement(transformation(extent={{-10,-60},{10,-40}})));
        SubstationTeststand.Components.Substation.Substation substation(
          m1_flow_nominal=m1_flow_nominal,
          m2_flow_nominal=m2_flow_nominal,
          dp1_nominal=dp1_nominal,
          dp2_nominal=dp2_nominal,
          dpValve_nominal=dpValve_nominal,
          hex(show_T=true))
          annotation (Placement(transformation(extent={{-10,20},{10,0}})));
        Components.SingleConsumerWithBooster
                                  singleConsumerWithBooster(
          m_flow_nominal=m2_flow_nominal,
          dp_nominal=dp2_nominal,
          T_supply_nominal=T_supply_secondary_nominal,
          T_return_nominal=T_return_secondary_nominal,
          Q_flow_nominal=Q_flow_nominal,
          scaling=2,
          fan(show_T=true, addPowerToMedium=false))
          annotation (Placement(transformation(extent={{-10,40},{10,60}})));

        Modelica.Blocks.Logical.Switch switch
          annotation (Placement(transformation(extent={{-44,6},{-36,14}})));
        Modelica.Blocks.Sources.Constant T_supply_secondary_high(k=
              T_supply_secondary_nominal) annotation (Placement(transformation(
              extent={{-4,4},{4,-4}},
              rotation=90,
              origin={-50,-4})));
        Modelica.Blocks.Sources.Constant T_supply_secondary_low(k=
              T_supply_secondary_nominal - 5.0) annotation (Placement(transformation(
              extent={{-4,-4},{4,4}},
              rotation=-90,
              origin={-50,24})));
      equation


        connect(substation.port_b2, singleConsumerWithBooster.port_a)
          annotation (Line(points={{-10,16},{-30,16},{-30,50},{-10,50}}, color=
                {0,127,255}));
        connect(singleConsumerWithBooster.port_b, substation.port_a2)
          annotation (Line(points={{10,50},{30,50},{30,16},{10,16}}, color={0,
                127,255}));
        connect(substation.port_b1, dHNetworkCtrl.port_a) annotation (Line(points={{10,
                4},{30,4},{30,-30},{5,-30},{5,-40}}, color={0,127,255}));
        connect(substation.port_a1, dHNetworkCtrl.port_b) annotation (Line(points={{-10,
                4},{-30,4},{-30,-30},{-5,-30},{-5,-40}}, color={0,127,255}));
        connect(T_supply_secondary_low.y, switch.u1) annotation (Line(points={{-50,19.6},
                {-50,13.2},{-44.8,13.2}}, color={0,0,127}));
        connect(T_supply_secondary_high.y, switch.u3) annotation (Line(points={{-50,0.4},
                {-50,6.8},{-44.8,6.8}}, color={0,0,127}));
        connect(switch.y, substation.T_supply_secondary_set)
          annotation (Line(points={{-35.6,10},{-12,10}}, color={0,0,127}));
        connect(switch.u2, dHNetworkCtrl.T_supply_secondary_ctrl) annotation (Line(
              points={{-44.8,10},{-70,10},{-70,-70},{0,-70},{0,-61}}, color={255,0,255}));
        annotation (
          Icon(
            coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
                  100}}), graphics={
              Ellipse(lineColor = {75,138,73},
                      fillColor={255,255,255},
                      fillPattern = FillPattern.Solid,
                      extent={{-100,-100},{100,100}}),
              Polygon(lineColor = {0,0,255},
                      fillColor = {75,138,73},
                      pattern = LinePattern.None,
                      fillPattern = FillPattern.Solid,
                      points={{-36,60},{64,0},{-36,-60},{-36,60}})}),
          Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                  {100,100}})),
          experiment(
            StopTime=172800,
            __Dymola_NumberOfIntervals=2880,
            __Dymola_Algorithm="Dassl"),
          Documentation(info="<html>
    <p>This model implements the thermal system of the JRA 1.1 multi-energy benchmark. The heating network is operated by controlling the pressure, keeping the pressure drop at the consumers above a certain threshold.</p>
    </html>"));
      end TestCase2_ReferenceSystem_subsystems;

      model TestCase2_Testbed

        parameter Modelica.SIunits.MassFlowRate m1_flow_nominal = 5 "Nominal mass flow rate primary side"
          annotation(Dialog(group="Nominal condition"));
        parameter Modelica.SIunits.MassFlowRate m2_flow_nominal = 5 "Nominal mass flow rate secondary side"
          annotation(Dialog(group="Nominal condition"));
        parameter Modelica.SIunits.PressureDifference dp1_nominal(displayUnit=
              "bar") = 200000 "Nominal pressure difference primary side"
              annotation(Dialog(group="Nominal condition"));
        parameter Modelica.SIunits.PressureDifference dp2_nominal(displayUnit=
              "bar") = 200000 "Nominal pressure difference secondary side"
              annotation(Dialog(group="Nominal condition"));
        parameter Modelica.SIunits.PressureDifference dpValve_nominal=200000
          "Nominal pressure drop of fully open valve" annotation(Dialog(group="Nominal condition"));

        parameter Modelica.SIunits.PressureDifference delta_p_secondary_nominal(displayUnit="bar") = 200000 "Nominal pressure difference for secondary side";
        parameter Modelica.SIunits.PressureDifference p_ref_secondary_nominal(displayUnit="bar") = IBPSA.Media.Water.p_default "Reference pressure level for secondary side";

        parameter Modelica.SIunits.Power Q_flow_nominal(displayUnit="kW") = 50000
          "Nominal heat flow rate" annotation(Dialog(group="Nominal condition"));
        parameter Modelica.SIunits.Temperature T_supply_primary_nominal=348.15
          "Nominal primary supply temperature" annotation(Dialog(group="Nominal condition"));
        parameter Modelica.SIunits.Temperature T_return_primary_nominal=308.15
          "Nominal primary return temperature" annotation(Dialog(group="Nominal condition"));
        parameter Modelica.SIunits.Temperature T_supply_secondary_nominal=341.15 "Nominal secondary supply temperature" annotation(Dialog(group="Nominal condition"));
        parameter Modelica.SIunits.Temperature T_return_secondary_nominal=301.15 "Nominal secondary return temperature" annotation(Dialog(group="Nominal condition"));

        Components.DHNetworkCtrl dHNetworkCtrl
          annotation (Placement(transformation(extent={{-10,-80},{10,-60}})));
        DisHeatLib.Utilities.ThermoHydraulicEquivalent.Demand_mT demand_mT(
            m_flow_nominal=4, dp_nominal=200000)
          annotation (Placement(transformation(extent={{-6,-46},{6,-34}})));
        DisHeatLib.Utilities.ThermoHydraulicEquivalent.Supply_pT supply_pT
          annotation (Placement(transformation(extent={{-6,24},{6,36}})));
        SubstationTeststand.TestStandSimplifiedPIDPressureCtrl
                                           testStandSimplifiedPIDPressureCtrl(
          m1_flow_nominal=m1_flow_nominal,
          m2_flow_nominal=m2_flow_nominal,
          dp1_nominal=dp1_nominal,
          dp2_nominal=delta_p_secondary_nominal,
          dpValve_nominal=dpValve_nominal,
          substation(hex(show_T=true)))
          annotation (Placement(transformation(extent={{-10,0},{10,-20}})));
        Modelica.Blocks.Sources.Constant delta_p_secondary(k=
              delta_p_secondary_nominal)
          annotation (Placement(transformation(extent={{80,-20},{60,0}})));
        Components.SingleConsumerWithBooster
                                  singleConsumerWithBooster(
          use_p_ref=false,
          m_flow_nominal=m2_flow_nominal,
          dp_nominal=delta_p_secondary_nominal,
          T_supply_nominal=T_supply_secondary_nominal,
          T_return_nominal=T_return_secondary_nominal,
          Q_flow_nominal=Q_flow_nominal,
          scaling=2,
          fan(show_T=true, addPowerToMedium=false))
          annotation (Placement(transformation(extent={{10,70},{-10,90}})));
        Modelica.Blocks.Logical.Switch switch
          annotation (Placement(transformation(extent={{-60,-20},{-40,0}})));
        Modelica.Blocks.Sources.Constant T_supply_secondary_high(k=
              T_supply_secondary_nominal) annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=90,
              origin={-70,-40})));
        Modelica.Blocks.Sources.Constant T_supply_secondary_low(k=
              T_supply_secondary_nominal - 5.0) annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=-90,
              origin={-70,20})));
      equation

        supply_pT.delta_p_set = 0;
        supply_pT.p_ref_set = p_ref_secondary_nominal;

        connect(demand_mT.port_b, dHNetworkCtrl.port_a) annotation (Line(points=
               {{6,-40},{26,-40},{26,-54},{5,-54},{5,-60}}, color={0,127,255}));
        connect(demand_mT.port_a, dHNetworkCtrl.port_b) annotation (Line(points=
               {{-6,-40},{-26,-40},{-26,-54},{-5,-54},{-5,-60}}, color={0,127,
                255}));
        connect(demand_mT.T_supply, testStandSimplifiedPIDPressureCtrl.T_supply_primary_set)
          annotation (Line(points={{-3.6,-46.6},{-3.6,-50},{-30,-50},{-30,-16},
                {-12,-16}}, color={0,0,127}));
        connect(supply_pT.T_return, testStandSimplifiedPIDPressureCtrl.T_return_secondary_set)
          annotation (Line(points={{-2.4,23.4},{-2.4,10},{-30,10},{-30,-4},{-12,
                -4}}, color={0,0,127}));
        connect(testStandSimplifiedPIDPressureCtrl.m_flow_return_secondary,
          supply_pT.m_flow) annotation (Line(points={{12,-10},{30,-10},{30,10},
                {2.4,10},{2.4,23.4}}, color={0,0,127}));
        connect(testStandSimplifiedPIDPressureCtrl.m_flow_return_primary,
          demand_mT.m_flow_set) annotation (Line(points={{-6,-21},{-6,-28},{-2.4,
                -28},{-2.4,-32.8}}, color={0,0,127}));
        connect(testStandSimplifiedPIDPressureCtrl.T_return_primary, demand_mT.T_return_set)
          annotation (Line(points={{6,-21},{6,-28},{2.4,-28},{2.4,-32.8}},
              color={0,0,127}));
        connect(testStandSimplifiedPIDPressureCtrl.T_supply_secondary,
          supply_pT.T_supply_set) annotation (Line(points={{6,1},{6,16},{12,16},
                {12,46},{3.6,46},{3.6,37.2}}, color={0,0,127}));
        connect(testStandSimplifiedPIDPressureCtrl.delta_p_secondary_set,
          delta_p_secondary.y) annotation (Line(points={{12,-4},{40,-4},{40,-10},
                {59,-10}}, color={0,0,127}));
        connect(testStandSimplifiedPIDPressureCtrl.delta_p_primary_set,
          demand_mT.delta_p) annotation (Line(points={{12,-16},{30,-16},{30,-50},
                {3.6,-50},{3.6,-46.6}}, color={0,0,127}));
        connect(singleConsumerWithBooster.port_a, supply_pT.port_b) annotation (
           Line(points={{10,80},{30,80},{30,30},{6,30}}, color={0,127,255}));
        connect(singleConsumerWithBooster.port_b, supply_pT.port_a) annotation (
           Line(points={{-10,80},{-32,80},{-32,30},{-6,30}}, color={0,127,255}));
        connect(T_supply_secondary_low.y, switch.u1) annotation (Line(points={{
                -70,9},{-70,-2},{-62,-2}}, color={0,0,127}));
        connect(T_supply_secondary_high.y, switch.u3) annotation (Line(points={
                {-70,-29},{-70,-18},{-62,-18}}, color={0,0,127}));
        connect(switch.y, testStandSimplifiedPIDPressureCtrl.T_supply_secondary_set)
          annotation (Line(points={{-39,-10},{-12,-10}}, color={0,0,127}));
        connect(dHNetworkCtrl.T_supply_secondary_ctrl, switch.u2) annotation (
            Line(points={{0,-81},{0,-90},{-88,-90},{-88,-10},{-62,-10}}, color=
                {255,0,255}));
        annotation (
          Icon(
            coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
                  100}}), graphics={
              Ellipse(lineColor = {75,138,73},
                      fillColor={255,255,255},
                      fillPattern = FillPattern.Solid,
                      extent={{-100,-100},{100,100}}),
              Polygon(lineColor = {0,0,255},
                      fillColor = {75,138,73},
                      pattern = LinePattern.None,
                      fillPattern = FillPattern.Solid,
                      points={{-36,60},{64,0},{-36,-60},{-36,60}})}),
          Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                  {100,100}})),
          experiment(
            StopTime=172800,
            __Dymola_NumberOfIntervals=2880,
            __Dymola_Algorithm="Dassl"),
          Documentation(info="<html>
    <p>This model implements the thermal system of the JRA 1.1 multi-energy benchmark. The heating network is operated by controlling the pressure, keeping the pressure drop at the consumers above a certain threshold.</p>
    </html>"));
      end TestCase2_Testbed;

      model TestCase2_Testbed_subsystems
        FMI.TeststandCtrl teststandCtrl
          annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
        FMI.SingleConsumerWithBooster consumer
          annotation (Placement(transformation(extent={{-10,50},{10,70}})));
        FMI.DHNetworkCtrl network
          annotation (Placement(transformation(extent={{-10,-70},{10,-50}})));
      equation
        connect(teststandCtrl.T_supply_secondary, consumer.T_supply_set)
          annotation (Line(points={{0,11},{0,20},{40,20},{40,90},{0,90},{0,72}},
              color={0,0,127}));
        connect(consumer.T_return, teststandCtrl.T_return_secondary_set)
          annotation (Line(points={{-6,49},{-6,40},{-40,40},{-40,6},{-12,6}},
              color={0,0,127}));
        connect(consumer.delta_p_secondary, teststandCtrl.delta_p_secondary_set)
          annotation (Line(points={{0,49},{0,34},{28,34},{28,5.8},{12,5.8}},
              color={0,0,127}));
        connect(consumer.m_flow, teststandCtrl.m_flow_return_secondary)
          annotation (Line(points={{6,49},{6,40},{34,40},{34,0},{12,0}}, color=
                {0,0,127}));
        connect(teststandCtrl.T_return_primary, network.T_return_set)
          annotation (Line(points={{4,-11},{4,-48}}, color={0,0,127}));
        connect(teststandCtrl.m_flow_return_primary, network.m_flow)
          annotation (Line(points={{-4,-11},{-4,-30},{-4,-48},{-3.8,-48}},
              color={0,0,127}));
        connect(network.delta_p, teststandCtrl.delta_p_primary_set) annotation (
           Line(points={{4,-71},{4,-80},{40,-80},{40,-6},{12,-6}}, color={0,0,
                127}));
        connect(network.T_supply, teststandCtrl.T_supply_primary_set)
          annotation (Line(points={{-3.8,-71},{-3.8,-80},{-34,-80},{-34,-6},{
                -12,-6}}, color={0,0,127}));
        connect(teststandCtrl.T_supply_secondary_ctrl, network.T_supply_secondary_ctrl)
          annotation (Line(points={{-12,0},{-40,0},{-40,-86},{0,-86},{0,-71}},
              color={255,0,255}));
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
          experiment(
            StopTime=172800,
            __Dymola_NumberOfIntervals=2880,
            __Dymola_Algorithm="Dassl"));
      end TestCase2_Testbed_subsystems;

    end TestCase2;
    annotation (
      Icon(
        graphics={
          Rectangle(
            lineColor={200,200,200},
            fillColor={248,248,248},
            fillPattern=FillPattern.HorizontalCylinder,
            extent={{-100,-100},{100,100}},
            radius=25),
          Polygon(
            origin={18,14},
            lineColor={78,138,73},
            fillColor={78,138,73},
            pattern=LinePattern.None,
            fillPattern=FillPattern.Solid,
            points={{-58.0,46.0},{42.0,-14.0},{-58.0,-74.0},{-58.0,46.0}}),
          Rectangle(
            extent={{-100,100},{100,-100}},
            lineColor={0,0,0},
            radius=25)}));
  end TestCases;

  package FMI
    model DHNetwork
      .TestbedExample.Components.DHNetwork dHNetwork(fileName=
            fileName, tableName=tableName)
        annotation (Placement(transformation(extent={{-10,-60},{10,-40}})));
      DisHeatLib.Utilities.ThermoHydraulicEquivalent.Demand_mT demand_mT(
          m_flow_nominal=4, dp_nominal=200000)
        annotation (Placement(transformation(extent={{-10,20},{10,40}})));
      Modelica.Blocks.Interfaces.RealInput m_flow annotation (Placement(
            transformation(
            extent={{-20,-20},{20,20}},
            rotation=-90,
            origin={-40,120}), iconTransformation(
            extent={{-20,-20},{20,20}},
            rotation=-90,
            origin={-38,120})));
      Modelica.Blocks.Interfaces.RealOutput T_supply annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=-90,
            origin={-40,-110}), iconTransformation(
            extent={{-10,-11},{10,11}},
            rotation=-90,
            origin={-38,-110})));
      Modelica.Blocks.Interfaces.RealInput T_return_set annotation (Placement(
            transformation(
            extent={{-20,-20},{20,20}},
            rotation=-90,
            origin={40,120}), iconTransformation(
            extent={{-20,-20},{20,20}},
            rotation=-90,
            origin={40,120})));
      Modelica.Blocks.Interfaces.RealOutput delta_p annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=-90,
            origin={40,-110}), iconTransformation(
            extent={{-10,-11},{10,11}},
            rotation=-90,
            origin={40,-110})));

      parameter String fileName="modelica://TestbedExample/data/network_data.txt"
        "File where matrix is stored" annotation(Dialog(group="Demand profile"),loadSelector(filter="Text files(*.txt);;CSV files (*.csv)",caption="Open data file"));

      parameter String tableName="NetworkData" "Table name on file" annotation(Dialog(group="Demand profile"));

    equation
      connect(demand_mT.port_b,dHNetwork. port_a) annotation (Line(points={{10,30},
              {20,30},{20,-20},{5,-20},{5,-40}},       color={0,127,255}));
      connect(m_flow, demand_mT.m_flow_set) annotation (Line(points={{-40,120},
              {-40,60},{-4,60},{-4,42}}, color={0,0,127}));
      connect(T_return_set, demand_mT.T_return_set) annotation (Line(points={{
              40,120},{40,60},{4,60},{4,42}}, color={0,0,127}));
      connect(demand_mT.delta_p, delta_p) annotation (Line(points={{6,19},{6,0},
              {40,0},{40,-110}}, color={0,0,127}));
      connect(demand_mT.port_a, dHNetwork.port_b) annotation (Line(points={{-10,
              30},{-20,30},{-20,-20},{-5,-20},{-5,-40}}, color={0,127,255}));
      connect(demand_mT.T_supply, T_supply) annotation (Line(points={{-6,19},{
              -6,0},{-40,0},{-40,-110}}, color={0,0,127}));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
                                           Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={0,0,0},
              fillColor={227,221,16},
              fillPattern=FillPattern.Solid,
              radius=20),
            Text(
              extent={{-60,-118},{60,-178}},
              lineColor={0,0,0},
              pattern=LinePattern.None,
              lineThickness=1,
              fillColor={0,0,173},
              fillPattern=FillPattern.Solid,
              textString="%name")}),                                 Diagram(
            coordinateSystem(preserveAspectRatio=false)));
    end DHNetwork;

    model DHNetworkCtrl
      Components.DHNetworkCtrl dHNetworkCtrl(fileName=fileName, tableName=
            tableName)
        annotation (Placement(transformation(extent={{-10,-60},{10,-40}})));
      DisHeatLib.Utilities.ThermoHydraulicEquivalent.Demand_mT demand_mT(
          m_flow_nominal=4, dp_nominal=200000)
        annotation (Placement(transformation(extent={{-10,20},{10,40}})));
      Modelica.Blocks.Interfaces.RealInput m_flow annotation (Placement(
            transformation(
            extent={{-20,-20},{20,20}},
            rotation=-90,
            origin={-40,120}), iconTransformation(
            extent={{-20,-20},{20,20}},
            rotation=-90,
            origin={-38,120})));
      Modelica.Blocks.Interfaces.RealOutput T_supply annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=-90,
            origin={-60,-110}), iconTransformation(
            extent={{-10,-11},{10,11}},
            rotation=-90,
            origin={-38,-110})));
      Modelica.Blocks.Interfaces.RealInput T_return_set annotation (Placement(
            transformation(
            extent={{-20,-20},{20,20}},
            rotation=-90,
            origin={40,120}), iconTransformation(
            extent={{-20,-20},{20,20}},
            rotation=-90,
            origin={40,120})));
      Modelica.Blocks.Interfaces.RealOutput delta_p annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=-90,
            origin={60,-110}), iconTransformation(
            extent={{-10,-11},{10,11}},
            rotation=-90,
            origin={40,-110})));
      Modelica.Blocks.Interfaces.BooleanOutput T_supply_secondary_ctrl
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=-90,
            origin={0,-110})));

      parameter String fileName="modelica://TestbedExample/data/network_data.txt"
        "File where matrix is stored" annotation(Dialog(group="Demand profile"),loadSelector(filter="Text files(*.txt);;CSV files (*.csv)",caption="Open data file"));

      parameter String tableName="NetworkData" "Table name on file" annotation(Dialog(group="Demand profile"));

    equation
      connect(demand_mT.port_b, dHNetworkCtrl.port_a) annotation (Line(points={
              {10,30},{20,30},{20,-20},{5,-20},{5,-40}}, color={0,127,255}));
      connect(m_flow, demand_mT.m_flow_set) annotation (Line(points={{-40,120},
              {-40,60},{-4,60},{-4,42}}, color={0,0,127}));
      connect(T_return_set, demand_mT.T_return_set) annotation (Line(points={{
              40,120},{40,60},{4,60},{4,42}}, color={0,0,127}));
      connect(demand_mT.delta_p, delta_p) annotation (Line(points={{6,19},{6,0},
              {60,0},{60,-110}}, color={0,0,127}));
      connect(demand_mT.port_a, dHNetworkCtrl.port_b) annotation (Line(points={
              {-10,30},{-20,30},{-20,-20},{-5,-20},{-5,-40}}, color={0,127,255}));
      connect(demand_mT.T_supply, T_supply) annotation (Line(points={{-6,19},{
              -6,0},{-60,0},{-60,-110}}, color={0,0,127}));
      connect(dHNetworkCtrl.T_supply_secondary_ctrl, T_supply_secondary_ctrl)
        annotation (Line(points={{0,-61},{0,-110}}, color={255,0,255}));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
                                           Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={0,0,0},
              fillColor={227,221,16},
              fillPattern=FillPattern.Solid,
              radius=20),
            Text(
              extent={{-60,-118},{60,-178}},
              lineColor={0,0,0},
              pattern=LinePattern.None,
              lineThickness=1,
              fillColor={0,0,173},
              fillPattern=FillPattern.Solid,
              textString="%name")}),                                 Diagram(
            coordinateSystem(preserveAspectRatio=false)));
    end DHNetworkCtrl;

    model Teststand

      parameter Modelica.SIunits.Temperature T_supply_primary_nominal=348.15 "Nominal supply temperature";
      parameter Modelica.SIunits.Temperature T_return_primary_nominal=308.15 "Nominal return temperature";

      parameter Modelica.SIunits.Temperature T_supply_secondary_nominal=341.15 "Setpoint for secondary supply temperature";

      parameter Modelica.SIunits.PressureDifference delta_p_secondary_nominal(displayUnit="bar") = 200000 "Nominal pressure difference for secondary side";

      Modelica.Blocks.Sources.Constant T_supply_secondary_set(k=
            T_supply_secondary_nominal)
        annotation (Placement(transformation(extent={{-60,-10},{-40,10}})));
      SubstationTeststand.TestStandSimplifiedPIDPressureCtrl testStandDummy(
        m1_flow_nominal=5,
        m2_flow_nominal=5,
        dp1_nominal=200000,
        dp2_nominal=delta_p_secondary_nominal,
        dpValve_nominal=400000,
        substation(hex(show_T=true)))
        annotation (Placement(transformation(extent={{-10,10},{10,-10}})));
      Modelica.Blocks.Interfaces.RealOutput m_flow_return_primary annotation (
          Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=-90,
            origin={-40,-110}), iconTransformation(
            extent={{-10,-11},{10,11}},
            rotation=-90,
            origin={-40,-110})));
      Modelica.Blocks.Interfaces.RealOutput T_return_primary annotation (
          Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=-90,
            origin={40,-110}), iconTransformation(
            extent={{-10,-11},{10,11}},
            rotation=-90,
            origin={40,-110})));
      Modelica.Blocks.Interfaces.RealInput T_return_secondary_set annotation (
          Placement(transformation(
            extent={{-20,-20},{20,20}},
            rotation=0,
            origin={-120,60}), iconTransformation(
            extent={{20,-20},{-20,20}},
            rotation=180,
            origin={-120,60})));
      Modelica.Blocks.Interfaces.RealInput T_supply_primary_set annotation (
          Placement(transformation(
            extent={{-20,-20},{20,20}},
            rotation=0,
            origin={-120,-60}), iconTransformation(
            extent={{20,-20},{-20,20}},
            rotation=180,
            origin={-120,-60})));
      Modelica.Blocks.Interfaces.RealInput m_flow_return_secondary annotation (
          Placement(transformation(
            extent={{20,-20},{-20,20}},
            rotation=0,
            origin={120,0}), iconTransformation(
            extent={{20,-20},{-20,20}},
            rotation=0,
            origin={120,0})));
      Modelica.Blocks.Interfaces.RealInput delta_p_secondary_set annotation (
          Placement(transformation(
            extent={{20,-20},{-20,20}},
            rotation=0,
            origin={120,60}), iconTransformation(
            extent={{20,-20},{-20,20}},
            rotation=0,
            origin={120,58})));
      Modelica.Blocks.Interfaces.RealInput delta_p_primary_set annotation (
          Placement(transformation(
            extent={{20,-20},{-20,20}},
            rotation=0,
            origin={120,-60}), iconTransformation(
            extent={{20,-20},{-20,20}},
            rotation=0,
            origin={120,-60})));
      Modelica.Blocks.Interfaces.RealOutput T_supply_secondary annotation (
          Placement(transformation(
            extent={{10,-10},{-10,10}},
            rotation=-90,
            origin={0,110}), iconTransformation(
            extent={{10,-11},{-10,11}},
            rotation=-90,
            origin={0,110})));
    equation
      connect(T_supply_secondary_set.y, testStandDummy.T_supply_secondary_set)
        annotation (Line(points={{-39,0},{-12,0}}, color={0,0,127}));
      connect(T_return_secondary_set, testStandDummy.T_return_secondary_set)
        annotation (Line(points={{-120,60},{-30,60},{-30,6},{-12,6}}, color={0,
              0,127}));
      connect(T_supply_primary_set, testStandDummy.T_supply_primary_set)
        annotation (Line(points={{-120,-60},{-28,-60},{-28,-6},{-12,-6}}, color=
             {0,0,127}));
      connect(delta_p_secondary_set, testStandDummy.delta_p_secondary_set)
        annotation (Line(points={{120,60},{30,60},{30,6},{12,6}}, color={0,0,
              127}));
      connect(m_flow_return_secondary, testStandDummy.m_flow_return_secondary)
        annotation (Line(points={{120,0},{12,0}}, color={0,0,127}));
      connect(delta_p_primary_set, testStandDummy.delta_p_primary_set)
        annotation (Line(points={{120,-60},{30,-60},{30,-6},{12,-6}}, color={0,
              0,127}));
      connect(m_flow_return_primary, testStandDummy.m_flow_return_primary)
        annotation (Line(points={{-40,-110},{-40,-80},{-6,-80},{-6,-11}}, color=
             {0,0,127}));
      connect(T_return_primary, testStandDummy.T_return_primary) annotation (
          Line(points={{40,-110},{40,-80},{6,-80},{6,-11}}, color={0,0,127}));
      connect(T_supply_secondary, testStandDummy.T_supply_secondary)
        annotation (Line(points={{0,110},{0,80},{6,80},{6,11}}, color={0,0,127}));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
                                           Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={0,0,0},
              fillColor={227,221,16},
              fillPattern=FillPattern.Solid,
              radius=20),
            Text(
              extent={{-60,-118},{60,-178}},
              lineColor={0,0,0},
              pattern=LinePattern.None,
              lineThickness=1,
              fillColor={0,0,173},
              fillPattern=FillPattern.Solid,
              textString="%name")}),                                 Diagram(
            coordinateSystem(preserveAspectRatio=false)));
    end Teststand;

    model TeststandCtrl

      parameter Modelica.SIunits.Temperature T_supply_primary_nominal=348.15 "Nominal supply temperature";
      parameter Modelica.SIunits.Temperature T_return_primary_nominal=308.15 "Nominal return temperature";

      parameter Modelica.SIunits.Temperature T_supply_secondary_nominal=341.15 "Setpoint for secondary supply temperature";

      parameter Modelica.SIunits.PressureDifference delta_p_secondary_nominal(displayUnit="bar") = 200000 "Nominal pressure difference for secondary side";

      SubstationTeststand.TestStandSimplifiedPIDPressureCtrl testStandDummy(
        m1_flow_nominal=5,
        m2_flow_nominal=5,
        dp1_nominal=200000,
        dp2_nominal=delta_p_secondary_nominal,
        dpValve_nominal=400000,
        substation(hex(show_T=true)))
        annotation (Placement(transformation(extent={{-10,10},{10,-10}})));
      Modelica.Blocks.Interfaces.RealOutput m_flow_return_primary annotation (
          Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=-90,
            origin={-40,-110}), iconTransformation(
            extent={{-10,-11},{10,11}},
            rotation=-90,
            origin={-40,-110})));
      Modelica.Blocks.Interfaces.RealOutput T_return_primary annotation (
          Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=-90,
            origin={40,-110}), iconTransformation(
            extent={{-10,-11},{10,11}},
            rotation=-90,
            origin={40,-110})));
      Modelica.Blocks.Interfaces.RealInput T_return_secondary_set annotation (
          Placement(transformation(
            extent={{-20,-20},{20,20}},
            rotation=0,
            origin={-120,60}), iconTransformation(
            extent={{20,-20},{-20,20}},
            rotation=180,
            origin={-120,60})));
      Modelica.Blocks.Interfaces.RealInput T_supply_primary_set annotation (
          Placement(transformation(
            extent={{-20,-20},{20,20}},
            rotation=0,
            origin={-120,-60}), iconTransformation(
            extent={{20,-20},{-20,20}},
            rotation=180,
            origin={-120,-60})));
      Modelica.Blocks.Interfaces.RealInput m_flow_return_secondary annotation (
          Placement(transformation(
            extent={{20,-20},{-20,20}},
            rotation=0,
            origin={120,0}), iconTransformation(
            extent={{20,-20},{-20,20}},
            rotation=0,
            origin={120,0})));
      Modelica.Blocks.Interfaces.RealInput delta_p_secondary_set annotation (
          Placement(transformation(
            extent={{20,-20},{-20,20}},
            rotation=0,
            origin={120,60}), iconTransformation(
            extent={{20,-20},{-20,20}},
            rotation=0,
            origin={120,58})));
      Modelica.Blocks.Interfaces.RealInput delta_p_primary_set annotation (
          Placement(transformation(
            extent={{20,-20},{-20,20}},
            rotation=0,
            origin={120,-60}), iconTransformation(
            extent={{20,-20},{-20,20}},
            rotation=0,
            origin={120,-60})));
      Modelica.Blocks.Interfaces.RealOutput T_supply_secondary annotation (
          Placement(transformation(
            extent={{10,-10},{-10,10}},
            rotation=-90,
            origin={0,110}), iconTransformation(
            extent={{10,-11},{-10,11}},
            rotation=-90,
            origin={0,110})));
      Modelica.Blocks.Interfaces.BooleanInput T_supply_secondary_ctrl
        annotation (Placement(transformation(extent={{-140,-20},{-100,20}})));
      Modelica.Blocks.Logical.Switch switch
        annotation (Placement(transformation(extent={{-60,-10},{-40,10}})));
      Modelica.Blocks.Sources.Constant T_supply_secondary_high(k=
            T_supply_secondary_nominal) annotation (Placement(transformation(
            extent={{-10,10},{10,-10}},
            rotation=90,
            origin={-70,-30})));
      Modelica.Blocks.Sources.Constant T_supply_secondary_low(k=
            T_supply_secondary_nominal - 5.0) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=-90,
            origin={-70,30})));
    equation
      connect(T_return_secondary_set, testStandDummy.T_return_secondary_set)
        annotation (Line(points={{-120,60},{-30,60},{-30,6},{-12,6}}, color={0,
              0,127}));
      connect(T_supply_primary_set, testStandDummy.T_supply_primary_set)
        annotation (Line(points={{-120,-60},{-28,-60},{-28,-6},{-12,-6}}, color=
             {0,0,127}));
      connect(delta_p_secondary_set, testStandDummy.delta_p_secondary_set)
        annotation (Line(points={{120,60},{30,60},{30,6},{12,6}}, color={0,0,
              127}));
      connect(m_flow_return_secondary, testStandDummy.m_flow_return_secondary)
        annotation (Line(points={{120,0},{12,0}}, color={0,0,127}));
      connect(delta_p_primary_set, testStandDummy.delta_p_primary_set)
        annotation (Line(points={{120,-60},{30,-60},{30,-6},{12,-6}}, color={0,
              0,127}));
      connect(m_flow_return_primary, testStandDummy.m_flow_return_primary)
        annotation (Line(points={{-40,-110},{-40,-80},{-6,-80},{-6,-11}}, color=
             {0,0,127}));
      connect(T_return_primary, testStandDummy.T_return_primary) annotation (
          Line(points={{40,-110},{40,-80},{6,-80},{6,-11}}, color={0,0,127}));
      connect(T_supply_secondary, testStandDummy.T_supply_secondary)
        annotation (Line(points={{0,110},{0,80},{6,80},{6,11}}, color={0,0,127}));
      connect(T_supply_secondary_low.y,switch. u1) annotation (Line(points={{-70,19},
              {-70,8},{-62,8}},         color={0,0,127}));
      connect(T_supply_secondary_high.y,switch. u3) annotation (Line(points={{-70,-19},
              {-70,-8},{-62,-8}},     color={0,0,127}));
      connect(switch.y, testStandDummy.T_supply_secondary_set)
        annotation (Line(points={{-39,0},{-12,0}}, color={0,0,127}));
      connect(T_supply_secondary_ctrl, switch.u2)
        annotation (Line(points={{-120,0},{-62,0}}, color={255,0,255}));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
                                           Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={0,0,0},
              fillColor={227,221,16},
              fillPattern=FillPattern.Solid,
              radius=20),
            Text(
              extent={{-60,-118},{60,-178}},
              lineColor={0,0,0},
              pattern=LinePattern.None,
              lineThickness=1,
              fillColor={0,0,173},
              fillPattern=FillPattern.Solid,
              textString="%name")}),                                 Diagram(
            coordinateSystem(preserveAspectRatio=false)));
    end TeststandCtrl;

    model SingleConsumer

      parameter Modelica.SIunits.MassFlowRate m_flow_nominal = 5 "Nominal mass flow rate secondary side"
        annotation(Dialog(group="Nominal condition"));
      parameter Modelica.SIunits.PressureDifference dp_nominal(displayUnit=
            "bar") = 200000 "Nominal pressure difference secondary side"
            annotation(Dialog(group="Nominal condition"));

      parameter Modelica.SIunits.Temperature T_supply_secondary_nominal=341.15 "Setpoint for secondary supply temperature";
      parameter Modelica.SIunits.Temperature T_return_secondary_nominal=301.15 "Setpoint for secondary supply temperature";

      parameter Modelica.SIunits.PressureDifference delta_p_secondary_nominal(displayUnit="bar") = 200000 "Nominal pressure difference for secondary side";
      parameter Modelica.SIunits.PressureDifference p_ref_secondary_nominal(displayUnit="bar") = IBPSA.Media.Water.p_default "Reference pressure level for secondary side";

      parameter Modelica.SIunits.Power Q_flow_nominal(displayUnit="kW") = 50000 "Nominal heat flow rate";

      parameter String fileName="modelica://TestbedExample/data/space_heating_load_profile_single_building_1week.txt"
        "File where matrix is stored" annotation(Dialog(group="Demand profile"),loadSelector(filter="Text files(*.txt);;CSV files (*.csv)",caption="Open data file"));
      parameter String tableName="HeatDemand" "Table name on file"  annotation(Dialog(group="Demand profile"));
      parameter Real scaling=2 "Scaling factor for heat demand" annotation(Dialog(group="Demand profile"));

      Modelica.Blocks.Sources.Constant delta_p_secondary_set(k=
            delta_p_secondary_nominal)
        annotation (Placement(transformation(extent={{-10,-10},{10,10}},
            rotation=-90,
            origin={0,-70})));
      DisHeatLib.Utilities.ThermoHydraulicEquivalent.Supply_pT supply_pT
        annotation (Placement(transformation(extent={{-10,-20},{10,0}})));
      Modelica.Blocks.Interfaces.RealOutput T_return annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=-90,
            origin={-60,-110}),
                              iconTransformation(
            extent={{-10,-11},{10,11}},
            rotation=-90,
            origin={-60,-110})));
      Modelica.Blocks.Interfaces.RealInput T_supply_set annotation (Placement(
            transformation(
            extent={{-20,-20},{20,20}},
            rotation=-90,
            origin={0,120}), iconTransformation(
            extent={{-20,-20},{20,20}},
            rotation=-90,
            origin={0,120})));
      Modelica.Blocks.Interfaces.RealOutput m_flow annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=-90,
            origin={60,-110}), iconTransformation(
            extent={{-10,-11},{10,11}},
            rotation=-90,
            origin={60,-110})));
      Modelica.Blocks.Interfaces.RealOutput delta_p_secondary annotation (
          Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=-90,
            origin={0,-110}),   iconTransformation(
            extent={{-10,-11},{10,11}},
            rotation=-90,
            origin={0,-110})));
      Components.SingleConsumer singleConsumer(
        use_p_ref=false,
        m_flow_nominal=m_flow_nominal,
        dp_nominal=delta_p_secondary_nominal,
        T_supply_nominal=T_supply_secondary_nominal,
        T_return_nominal=T_return_secondary_nominal,
        Q_flow_nominal=Q_flow_nominal,
        fileName=fileName,
        tableName=tableName,
        scaling=scaling,
        fan(show_T=true, addPowerToMedium=false))
        annotation (Placement(transformation(extent={{10,40},{-10,60}})));
    equation

      supply_pT.delta_p_set = 0;
      supply_pT.p_ref_set = p_ref_secondary_nominal;

      connect(T_supply_set, supply_pT.T_supply_set) annotation (Line(points={{0,120},
              {0,80},{20,80},{20,20},{6,20},{6,2}},          color={0,0,127}));
      connect(T_return, supply_pT.T_return) annotation (Line(points={{-60,-110},{-60,
              -40},{-4,-40},{-4,-21}},   color={0,0,127}));
      connect(m_flow, supply_pT.m_flow) annotation (Line(points={{60,-110},{60,-40},
              {4,-40},{4,-21}},      color={0,0,127}));
      connect(delta_p_secondary, delta_p_secondary)
        annotation (Line(points={{0,-110},{0,-110}},     color={0,0,127}));
      connect(delta_p_secondary, delta_p_secondary_set.y) annotation (Line(
            points={{0,-110},{0,-81},{-1.9984e-15,-81}},
                                                 color={0,0,127}));
      connect(singleConsumer.port_a, supply_pT.port_b) annotation (Line(points={{10,
              50},{30,50},{30,-10},{10,-10}}, color={0,127,255}));
      connect(supply_pT.port_a, singleConsumer.port_b) annotation (Line(points={{-10,
              -10},{-28,-10},{-28,50},{-10,50}}, color={0,127,255}));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
                                           Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={0,0,0},
              fillColor={227,221,16},
              fillPattern=FillPattern.Solid,
              radius=20),
            Text(
              extent={{-60,-118},{60,-178}},
              lineColor={0,0,0},
              pattern=LinePattern.None,
              lineThickness=1,
              fillColor={0,0,173},
              fillPattern=FillPattern.Solid,
              textString="%name")}),                                 Diagram(
            coordinateSystem(preserveAspectRatio=false)));
    end SingleConsumer;

    model SingleConsumerWithBooster

      parameter Modelica.SIunits.MassFlowRate m_flow_nominal = 5 "Nominal mass flow rate secondary side"
        annotation(Dialog(group="Nominal condition"));
      parameter Modelica.SIunits.PressureDifference dp_nominal(displayUnit=
            "bar") = 200000 "Nominal pressure difference secondary side"
            annotation(Dialog(group="Nominal condition"));

      parameter Modelica.SIunits.Temperature T_supply_secondary_nominal=341.15 "Setpoint for secondary supply temperature";
      parameter Modelica.SIunits.Temperature T_return_secondary_nominal=301.15 "Setpoint for secondary supply temperature";

      parameter Modelica.SIunits.PressureDifference delta_p_secondary_nominal(displayUnit="bar") = 200000 "Nominal pressure difference for secondary side";
      parameter Modelica.SIunits.PressureDifference p_ref_secondary_nominal(displayUnit="bar") = IBPSA.Media.Water.p_default "Reference pressure level for secondary side";

      parameter Modelica.SIunits.Power Q_flow_nominal(displayUnit="kW") = 50000 "Nominal heat flow rate";

      parameter String fileName="modelica://TestbedExample/data/space_heating_load_profile_single_building_1week.txt"
        "File where matrix is stored" annotation(Dialog(group="Demand profile"),loadSelector(filter="Text files(*.txt);;CSV files (*.csv)",caption="Open data file"));
      parameter String tableName="HeatDemand" "Table name on file"  annotation(Dialog(group="Demand profile"));
      parameter Real scaling=2 "Scaling factor for heat demand" annotation(Dialog(group="Demand profile"));

      Modelica.Blocks.Sources.Constant delta_p_secondary_set(k=
            delta_p_secondary_nominal)
        annotation (Placement(transformation(extent={{-10,-10},{10,10}},
            rotation=-90,
            origin={0,-70})));
      DisHeatLib.Utilities.ThermoHydraulicEquivalent.Supply_pT supply_pT
        annotation (Placement(transformation(extent={{-10,-20},{10,0}})));
      Modelica.Blocks.Interfaces.RealOutput T_return annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=-90,
            origin={-60,-110}),
                              iconTransformation(
            extent={{-10,-11},{10,11}},
            rotation=-90,
            origin={-60,-110})));
      Modelica.Blocks.Interfaces.RealInput T_supply_set annotation (Placement(
            transformation(
            extent={{-20,-20},{20,20}},
            rotation=-90,
            origin={0,120}), iconTransformation(
            extent={{-20,-20},{20,20}},
            rotation=-90,
            origin={0,120})));
      Modelica.Blocks.Interfaces.RealOutput m_flow annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=-90,
            origin={60,-110}), iconTransformation(
            extent={{-10,-11},{10,11}},
            rotation=-90,
            origin={60,-110})));
      Modelica.Blocks.Interfaces.RealOutput delta_p_secondary annotation (
          Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=-90,
            origin={0,-110}),   iconTransformation(
            extent={{-10,-11},{10,11}},
            rotation=-90,
            origin={0,-110})));
      Components.SingleConsumerWithBooster
                                singleConsumerWithBooster(
        use_p_ref=false,
        m_flow_nominal=m_flow_nominal,
        dp_nominal=delta_p_secondary_nominal,
        T_supply_nominal=T_supply_secondary_nominal,
        T_return_nominal=T_return_secondary_nominal,
        Q_flow_nominal=Q_flow_nominal,
        fileName=fileName,
        tableName=tableName,
        scaling=scaling,
        fan(show_T=true, addPowerToMedium=false))
        annotation (Placement(transformation(extent={{10,40},{-10,60}})));
    equation

      supply_pT.delta_p_set = 0;
      supply_pT.p_ref_set = p_ref_secondary_nominal;

      connect(T_supply_set, supply_pT.T_supply_set) annotation (Line(points={{0,120},
              {0,80},{20,80},{20,20},{6,20},{6,2}},          color={0,0,127}));
      connect(T_return, supply_pT.T_return) annotation (Line(points={{-60,-110},{-60,
              -40},{-4,-40},{-4,-21}},   color={0,0,127}));
      connect(m_flow, supply_pT.m_flow) annotation (Line(points={{60,-110},{60,-40},
              {4,-40},{4,-21}},      color={0,0,127}));
      connect(delta_p_secondary, delta_p_secondary)
        annotation (Line(points={{0,-110},{0,-110}},     color={0,0,127}));
      connect(delta_p_secondary, delta_p_secondary_set.y) annotation (Line(
            points={{0,-110},{0,-81},{-1.9984e-15,-81}},
                                                 color={0,0,127}));
      connect(singleConsumerWithBooster.port_a, supply_pT.port_b) annotation (
          Line(points={{10,50},{30,50},{30,-10},{10,-10}}, color={0,127,255}));
      connect(supply_pT.port_a, singleConsumerWithBooster.port_b) annotation (
          Line(points={{-10,-10},{-28,-10},{-28,50},{-10,50}}, color={0,127,255}));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
                                           Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={0,0,0},
              fillColor={227,221,16},
              fillPattern=FillPattern.Solid,
              radius=20),
            Text(
              extent={{-60,-118},{60,-178}},
              lineColor={0,0,0},
              pattern=LinePattern.None,
              lineThickness=1,
              fillColor={0,0,173},
              fillPattern=FillPattern.Solid,
              textString="%name")}),                                 Diagram(
            coordinateSystem(preserveAspectRatio=false)));
    end SingleConsumerWithBooster;
    annotation (Icon(graphics={
          Rectangle(
            lineColor={200,200,200},
            fillColor={248,248,248},
            fillPattern=FillPattern.HorizontalCylinder,
            extent={{-100,-100},{100,100}},
            radius=25),
          Rectangle(
            extent={{-100,100},{100,-100}},
            lineColor={0,0,0},
            radius=25),
          Bitmap(extent={{-64,-70},{86,64}}, fileName=
                "modelica://TestbedExample/images/fmi.svg")}));
  end FMI;
  annotation (uses(
      DisHeatLib(version="1.2"),
      IBPSA(version="3.0.0"),
      Modelica(version="3.2.3")));
end TestbedExample;
