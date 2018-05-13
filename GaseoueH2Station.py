
import CoolProp.CoolProp as CP
import numpy as np
import pandas as pd
import math
import json


def compressor_energy_demand(p_init, compressor_stage, prate_cascade, e_isentropic):
    """Compressor Energy Demand, Assumptions: T_in_2+=40°C"""

    w_spec_cascade = 0
    p_out_last_stage = p_init*100000
    for i in range(1, compressor_stage+1):
        if i == 1:
            T_in_stage = 298.15  # K
        else:
            T_in_stage = 313.15  # K

        p_in_stage = p_out_last_stage
        p_out_stage = p_in_stage*prate_cascade
        entropy_in = CP.PropsSI('S', 'T', T_in_stage,
                                'P', p_in_stage, "hydrogen")
        h_in = CP.PropsSI('H', 'S', entropy_in, 'P', p_in_stage, "hydrogen")
        h_out_is = CP.PropsSI('H', 'S', entropy_in, 'P',
                              p_out_stage, "hydrogen")
        h_out = (h_out_is-h_in)/e_isentropic+h_in
        w_comp = h_out - h_in
        p_out_last_stage = p_out_stage

        w_spec_cascade += w_comp/3600000

    return w_spec_cascade


def booster_compressor_energy_demand(p_init, prate_booster, e_isentropic):
    """Booster Compressor Energy Demand, Assumptions: T_in_2+=40°C"""

    w_spec_booster = 0
    p_out_last_stage = p_init*100000
    for i in range(1, 3):
        if i == 1:
            T_in_stage = 298.15  # K
        else:
            T_in_stage = 313.15  # K

        p_in_stage = p_out_last_stage
        p_out_stage = p_in_stage*prate_booster
        entropy_in = CP.PropsSI('S', 'T', T_in_stage,
                                'P', p_in_stage, "hydrogen")
        h_in = CP.PropsSI('H', 'S', entropy_in, 'P', p_in_stage, "hydrogen")
        h_out_is = CP.PropsSI('H', 'S', entropy_in, 'P',
                              p_out_stage, "hydrogen")
        h_out = (h_out_is-h_in)/e_isentropic+h_in
        w_comp = h_out - h_in
        p_out_last_stage = p_out_stage

        w_spec_booster += w_comp/3600000

    return w_spec_booster


class GaseousH2Station():

    """Design of gaseous h2 station"""

    def __init__(self, senario):

        # Station Type
        self.station_type = senario['station type']
        # Vehicle Service Pressure[bar]
        self.vehicle_service_pressure = senario['vehicle service pressure']
        # Refueling Station Size [kg/day]
        self.station_size = senario['station size']
        # Hydrogen Source
        self.hydrogen_source = senario['hydrogen source']
        # Dispensing Options to vehicle Tank
        self.dispensing_option = senario['dispensing option']
        # Compressor Stages
        self.compressor_stage = senario['compressor stage']
        # Vehicle Lingering time (min)
        self.t_linger = senario['vehicle Lingering time']
        # Utilization Rate
        self.utilization_rate = senario['utilization rate']

        # total vehicle fill time (min)
        self.t_fill = senario['total vehicle fill time']
        # Max. Dispensed Amount per Vehicle (kg)
        self.m_dispensed = senario['max. dispensed amount per vehicle']

        # Maximum number of vehicle back-to-back fills per HOSE during peak hours
        self.n_max = 60/(self.t_linger+self.t_fill)
        # Number of Hoses
        self.n_hoses = 1 if self.station_size == 212 else (
            2 if self.station_size == 420 else 4)

    def get_p_min(self):

        if self.hydrogen_source == "20 bar H2 supply Pipline":
            # else input("Please the Pressure of Hydrogen Delivered to Refueling Station (bar)")
            p_delivered = 20
            p_min = p_delivered

        elif self.hydrogen_source == "Tube-trailer supply":
            # input("Design Input Value:Please enter the maximum pressure in the Tube Trailer (bar): ") #Maximum Pressure in the Tube Trailer (bar)
            p_trailer_max = 200
            # input("Design Input Value:Please enter the minmum pressure in the Tube Trailer (bar): ") #Minimum Pressure in the Tube Trailer (bar)
            p_trailer_min = 50
            p_min = p_trailer_min

        else:
            print('hydrogen source is wrong!')
            p_min = 20
        return p_min

    def low_pressure_storage_design(self):
        p_storage_max = 250  # bar
        p_storage_min = 50  # bar
        # operating temperature of th station
        T_oper = 25  # ℃
        n_storage = 1
        Density_storage_max = CP.PropsSI(
            "D", "T", T_oper+273.15, "P", p_storage_max*100000, "hydrogen")
        Density_storage_min = CP.PropsSI(
            "D", "T", T_oper+273.15, "P", p_storage_min*100000, "hydrogen")
        percentage_useable = (Density_storage_max -
                              Density_storage_min)/Density_storage_max
        capacity_storage = self.station_size/percentage_useable
        low_pressure_storage_design = {'n_storage': n_storage,
                                       'capacity_storage': capacity_storage,
                                       'percentage_useable': percentage_useable}
        return low_pressure_storage_design

    def get_input_dispenser(self):
        """Dispenser Input"""
        # Vehicle Tank maximum pressure at end of dispensing [bar]#At 25oC
        p_max_dispenser = 1.25*self.vehicle_service_pressure

        if self.dispensing_option == "Cascade Dispensing":
            t_cascade = self.t_fill  # Vehicle fill time through cascade (min)
            t_booster = 0
        else:
            # Vehicle fill time through cascade (min)
            t_cascade = self.t_fill-2
            t_booster = 2  # Vehicle fill time through booster compressor (min)
        # self.t_linger #Vehicle Lingering time (min)
        # self.n_hoses #Hoses per dispenser
        # self.n_max #Maximum Dispensed Amount per Vehicle (kg)
        p_vessels = self.get_p_vessels()
        p_hp_switch = p_vessels.loc['switch', 'hp']
        v_tank_vehicle = self.n_hoses*self.m_dispensed / \
            CP.PropsSI("D", "P", p_hp_switch*100000, "T", 25 +
                       273.15, "hydrogen")  # Vehicle tank size (m³)
        input_dispenser = {'p_max_dispenser': p_max_dispenser,
                           't_cascade': t_cascade,
                           't_booster': t_booster,
                           'v_tank_vehicle': v_tank_vehicle}

        return input_dispenser

    def get_p_vessels(self):
        """ to get the pressure of Storage Vessels"""
        type = ['max', 'min', 'switch']
        # lp=low pressure, mp= middle pressure, hp=high pressure
        level = ['lp', 'mp', 'hp']

        p_vessels = pd.DataFrame(np.arange(9).reshape(
            (3, 3)), index=type, columns=level)  # create a dataframe

        # Maximum Pressure service pressure= Cascade Dispensing=vehicle service pressure*14,5*1,25+1000;350*14,5*1,25+1000)
        p_vessels.loc['max', :] = self.vehicle_service_pressure*1.25 + 70 \
            if self.dispensing_option == "Cascade Dispensing" else 350*1.25+70
        # for low pressure vessel: Minimum Pressure=0,35*maximum pressure
        p_vessels.loc['min', 'lp'] = 0.35*p_vessels.loc['max', 'lp']
        # for middle pressure vessel: Minimum Pressure=0,65*maximum pressure
        p_vessels.loc['min', 'mp'] = 0.65*p_vessels.loc['max', 'mp']
        # for hp pressure vessel: Minimum Pressure=0,85*maximum pressure
        p_vessels.loc['min', 'hp'] = 0.85*p_vessels.loc['max', 'hp']
        p_vessels.loc['switch', 'lp'] = 0.906 * \
            p_vessels.loc['min', 'lp']  # * 0,871;0,895; 0,906
        p_vessels.loc['switch', 'mp'] = 0.895*p_vessels.loc['min', 'mp']
        p_vessels.loc['switch', 'hp'] = 0.871*p_vessels.loc['min', 'hp']

        return p_vessels

    def refrigeration_design(self):
        """ Refrigeration Unit Energy demand [kWh/kg]"""
        if self.vehicle_service_pressure == 350:
            T_dispensing_max = -20
        else:
            T_dispensing_max = -40
        p_max_dispenser = self.get_input_dispenser()['p_max_dispenser']
        T_ref_design = 35
        T_amb = 25
        COP = 1.6*math.e**(-0.018*T_amb)
        h_1 = CP.PropsSI('H', 'T', T_ref_design+273.15,
                         'P', p_max_dispenser*100000, "hydrogen")
        h_2 = CP.PropsSI('H', 'T', T_dispensing_max+273.15,
                         'P', p_max_dispenser*100000, "hydrogen")
        # Refrigeration Unit Energy demand [kWh/kg]
        w_refuelling = (h_1 - h_2)/COP/3600000
        # HX mass
        m_HX = w_refuelling*3600*self.m_dispensed/5/0.837 ##5K allowed temperature rise
        # Power of refrigeration unit
        P_ref = w_refuelling*self.m_dispensed*self.n_max
        W_overhead = 25*math.log(T_amb)-21  # KWh/day
        w_overhead = W_overhead/self.station_size/self.utilization_rate  # KWh/kg
        w_refrigeration = w_refuelling+w_overhead
        n_ref = self.n_hoses

        refrigeration_design = {'n_ref': n_ref,
                                'T_dispensing_max': T_dispensing_max,
                                'w_refuelling': w_refuelling,
                                'w_overhead': w_overhead,
                                'w_refrigeration': w_refrigeration,
                                'm_HX': m_HX,
                                'P_ref': P_ref}
        return refrigeration_design

    def get_input_compressor(self):
        """Compressor Input"""
        # SAE J2601   #input("Compressor Input: Please enter the Compressor Discharge Temperature Before Intercooling (C): ")
        T_comp_before_cooling = 150
        # SAE J2601    #input("Compressor Input: Please enter Maximum vehicle tank temperature (C): ")
        T_vehicle_tank_max = 85
        # Operating Storage Temperature (degrees C)   #input("Compressor Input: Please enter Operating Storage Temperature (degrees C): ")
        T_oper = 25
        T_comp_after_cooling = 40
        # Main Compressor Discharge Pressure (bar),Corresponds to cascade vessels maximum pressure
        p_discharge = self.get_p_vessels().loc['max', 'hp']+25
        # self.compressor_stages
        # Isentropic Compressor Efficiency for Refueling Station Compressor (%)
        e_isentropic = 0.8
        f_sizing = 1.1  # Compressor Motor Sizing Factor, H2A base case value is 110% based on typical industrial practices
        # Hydrogen Lost During Compression (% of Feed H2)
        e_lost_compressor = 0.005
        input_compressor = {'p_discharge': p_discharge,
                            'e_isentropic': e_isentropic,
                            'f_sizing': f_sizing,
                            'e_lost_compressor': e_lost_compressor}

        return input_compressor

    def compressor_design(self):
        """ Compressor Design"""
        input_compressor = self.get_input_compressor()

        # Cm, Minimum Compressor Capacity (kg/hr)
        m_compressor_desired = self.n_hoses * \
            self.m_dispensed*60/(self.t_linger+self.t_fill)
        # Number of Compressors in Operation at Any Time
        n_compressor = math.ceil(m_compressor_desired/35)
        # Design Compressor Flow Rate (kg/hr)
        m_compressor_design = m_compressor_desired/n_compressor
        # n_dispensers_required = self.n_hoses
        p_min = self.get_p_min()  # get p_min
        prate_compressor = 10**((math.log10(input_compressor['p_discharge'])-math.log10(
            p_min))/self.compressor_stage)  # Pressure Ratio per stage cascade compressor
        # specific energy requirement  for the compression (kWh/kg)
        w_spec_compressor = compressor_energy_demand(
            p_min, self.compressor_stage, prate_compressor, input_compressor['e_isentropic'])
        # Theoretical Power Requirement for Each Compressor (kW)
        P_compressor_theo = w_spec_compressor*m_compressor_design
        # Actual Shaft Power Requirement for Each Compressor (kW)
        P_compressor_actual = P_compressor_theo
        e_motor = (0.00008*(math.log(P_compressor_actual))**4-0.0015*(math.log(P_compressor_actual))
                   ** 3+0.0061*(math.log(P_compressor_actual))**2+0.0311*(math.log(P_compressor_actual))+0.7617)
        P_compressor_motor = P_compressor_actual / \
            e_motor*input_compressor['f_sizing']
        # specific total energy requirement  (kWh/kg)
        w_sepc_compressor_motor = w_spec_compressor/e_motor

        compresor_design = {'m_compressor_design': m_compressor_design,
                            'prate_compressor': prate_compressor,
                            'n_compressor': n_compressor,
                            'P_compressor_actual': P_compressor_actual,
                            'P_compressor_motor': P_compressor_motor,
                            'w_sepc_compressor_motor': w_sepc_compressor_motor}
        return compresor_design

    def get_input_booster_compressor(self):
        """ Booster Compressor Input"""
        e_isentropic = 0.8
        # Number of Booster Stages
        booster_stage = 2 if self.dispensing_option == "Booster Compressor Dispensing" else 0
        # Booster Compressor Discharge Pressure (bar)
        p_booster_discharge = self.vehicle_service_pressure * \
            1.25 + \
            25 if self.dispensing_option == "Booster Compressor Dispensing" else 0

        input_booster_compressor = {'p_booster_discharge': p_booster_discharge,
                                    'e_isentropic': e_isentropic,
                                    'booster_stage': booster_stage}

        return input_booster_compressor

    def booster_compressor_design(self):
        """ Booster Compressor Design"""
        input_booster_compressor = self.get_input_booster_compressor()
        p_booster_discharge = input_booster_compressor['p_booster_discharge']
        e_isentropic = input_booster_compressor['e_isentropic']
        booster_stage = input_booster_compressor['booster_stage']
        # Booster compressor throughput [kg/min]
        m_booster = 42
        # Design Booster Compressor Flow Rate (kg/hr)
        m_booster_desired = self.m_dispensed*60 / \
            self.t_fill if self.dispensing_option == "Booster Compressor Dispensing" else 0
        n_booster = self.n_hoses*math.ceil(m_booster_desired / m_booster)
        m_booster_design = m_booster_desired*self.n_hoses / \
            n_booster if self.dispensing_option == "Booster Compressor Dispensing" else 0
        # Isentropic Compressor Efficiency for Refueling Station Compressor (%)
        p_hp_vessels_min = self.get_p_vessels().loc['min', 'hp']
        # Pressure Ratio per stage Booster compressor
        prate_booster = 10**((math.log10(p_booster_discharge) -
                              math.log10(p_hp_vessels_min))/booster_stage) if self.dispensing_option == "Booster Compressor Dispensing" else 0
        # Booster Compressor Shaft Specific Power (kWh/kg)
        w_spec_booster = booster_compressor_energy_demand(
            p_hp_vessels_min, prate_booster, e_isentropic) if self.dispensing_option == "Booster Compressor Dispensing" else 0
        P_shaft_booster = w_spec_booster * \
            m_booster_design  # Booster Shaft Power (kW)

        if P_shaft_booster != 0:
            e_motor_booster = (0.00008*(math.log(P_shaft_booster))**4-0.0015*(math.log(P_shaft_booster))**3
                               + 0.0061*(math.log(P_shaft_booster))**2+0.0311*(math.log(P_shaft_booster))+0.7617)
        else:
            e_motor_booster = 1  # Motor Efficiency for Booster Compressor
        # Theoretical Power Requirement for Booster Compressor (kW)
        P_booster_theo = P_shaft_booster
        # Motor Rating per Booster Compressor (kW)
        P_booster_motor = P_booster_theo / \
            e_motor_booster if self.dispensing_option == "Booster Compressor Dispensing" else 0
        # specific total energy requirement  (kWh/kg)
        w_sepc_booster_motor = w_spec_booster / \
            e_motor_booster if self.dispensing_option == "Booster Compressor Dispensing" else 0

        booster_compresor_design = {'n_booster': n_booster,
                                    'm_booster_design': m_booster_design,
                                    'prate_booster': prate_booster,
                                    'P_booster_motor': P_booster_motor,
                                    'w_sepc_booster_motor': w_sepc_booster_motor}
        return booster_compresor_design

    def cascade_calculation(self):
        """Design of Cascade System"""
        # t_oper = 24  #Hours the Refueling Station Is Operating(h)
        # m_hour_ava = self.m_hour_ava  #Average Hourly Demand During a Day (kg/hr)
        T_oper = 25  # ℃
        v_tank_vehicle = self.get_input_dispenser()['v_tank_vehicle']
        n_cascade_system = 1  # Optimum Number of Cascade Systems
        v_cascade_per_storage = 0.05  # Volume of Cascade storage [m³]
        variable = ['density_at_pmax', 'density_at_pmin', 'delta_density',  # list of needed variates
                    'm_vehicle_tank', 'm_cascade_vessel', 'v_min_cascade',
                    'n_tank', 'v_cascade', 'capacity_at_pmax',
                    'capacity_at_pmin', 'delta_capacity', 'percent_cascade_useable']
        # lp=low pressure mp= middle pressure hp=high pressure
        level = ['lp', 'mp', 'hp']
        p_vessels = self.get_p_vessels()
        cascade_design = pd.DataFrame(np.arange(36).reshape(
            (12, 3)), index=variable, columns=level)
        for i in range(3):
            # ### Low Pressure Cascade Storage Vessel
            # Density at Maximum Specified Unit Pressure (kg/m³)
            cascade_design.iloc[0, i] = CP.PropsSI(
                "D", "T", T_oper+273.15, "P", p_vessels.iloc[0, i]*100000, "hydrogen")
            # Density at Minimum Specified Unit Pressure (kg/m³)
            cascade_design.iloc[1, i] = CP.PropsSI(
                "D", "T", T_oper+273.15, "P", p_vessels.iloc[1, i]*100000, "hydrogen")
            # Difference Between Density at Maximum and Minimum Specified Unit Pressure - Cascade Tank (kg/m³)
            cascade_design.iloc[2, i] = cascade_design.iloc[0,
                                                            i] - cascade_design.iloc[1, i]
            # Maximum mass (kg) inside vehicle tank
            cascade_design.iloc[3, i] = v_tank_vehicle*CP.PropsSI(
                "D", "T", T_oper+273.15, "P", p_vessels.iloc[2, i]*100000, "hydrogen")
            # Mass required from  cascade vessel (kg)
            cascade_design.iloc[4, i] = cascade_design.iloc[3, i] if i == 0 else (
                cascade_design.iloc[3, i]-cascade_design.iloc[3, i-1])
            # Minimum Cascade Volume (m³)
            cascade_design.iloc[5, i] = cascade_design.iloc[4, i] / \
                cascade_design.iloc[2, i]
            # Number of Storage Tanks
            cascade_design.iloc[6, i] = math.ceil(
                cascade_design.iloc[5, i]/v_cascade_per_storage)
            # Cascade Volume [m³]
            cascade_design.iloc[7, i] = cascade_design.iloc[6,
                                                            i]*v_cascade_per_storage
            # Tank Capacity at Maximum Specified Unit Pressure (kg)
            cascade_design.iloc[8, i] = cascade_design.iloc[7, i]*CP.PropsSI(
                "D", "T", T_oper+273.15, "P", p_vessels.iloc[0, i]*100000, "hydrogen")
            # Tank Capacity at Minimum Specified Unit Pressure (kg)
            cascade_design.iloc[9, i] = cascade_design.iloc[7, i]*CP.PropsSI(
                "D", "T", T_oper+273.15, "P", p_vessels.iloc[1, i]*100000, "hydrogen")
            # Difference Between Capacity at Maximum and Minimum Specified Unit Pressure -  Cascade Tank (kg)
            cascade_design.iloc[10, i] = cascade_design.iloc[8,
                                                             i]-cascade_design.iloc[9, i]
            # Useable Capacity at Peak Demane
            cascade_design.iloc[11, i] = cascade_design.iloc[10,
                                                             i]/cascade_design.iloc[8, i]

        # Maximum Dispensible Amount from Cascade (kg)
        m_dispensiable = cascade_design.iloc[10, :].sum()
        # Number of Storage Tanks
        n_tank_total = cascade_design.iloc[6, :].sum()
        # Cascade Vessel Capacity/Unit (kg)
        capacity_at_pmax_sum = cascade_design.iloc[8, :].sum()/n_cascade_system
        # Cascade Size as a Percent of Average Daily Demand
        percent_cascad_of_demand = n_cascade_system*capacity_at_pmax_sum / \
            self.station_size
        cascade_system = {'cascade_design': cascade_design,
                          'm_dispensiable': m_dispensiable,
                          'n_tank_total': n_tank_total,
                          'capacity_at_pmax_sum': capacity_at_pmax_sum,
                          'percent_cascad_of_demand': percent_cascad_of_demand}
        return cascade_system


class Investment():
    """
    to calculat all the investment
    """
    def __init__(self,senario):
        
        self.station = GaseousH2Station(senario)
        self.f_install = 1.3
        self.rate_exchange = 0.82  # 2018 conversion of exchange rate

    def cost_lp_storage(self):
        """Investment of Low-Pressure Storage for Hourly Surge"""
        low_pressure_storage_design = self.station.low_pressure_storage_design()
        n_storage = low_pressure_storage_design['n_storage']
        capacity_storage = low_pressure_storage_design['capacity_storage']
        cost_lp_storage = n_storage*capacity_storage*1200*self.f_install*self.rate_exchange

        return cost_lp_storage

    def cost_dispenser(self):
        """"Investment of Dispenser"""
        cost_dispenser = self.station.n_hoses * 65000*self.f_install*self.rate_exchange

        return cost_dispenser

    def cost_refrigeration(self):
        """Investment of Refrigeration Equipment"""
        refrigeration_design = self.station.refrigeration_design()
        n_ref = refrigeration_design['n_ref']
        P_ref = refrigeration_design['P_ref']
        m_HX = refrigeration_design['m_HX']
        T_dispensing_max = refrigeration_design['T_dispensing_max']
        if m_HX > 1000:
            cost_refrigeration = n_ref*(14000*(28.43*P_ref/(T_dispensing_max+273.15))
                                        ** 0.8579+35000*(m_HX/1000)**0.9)*self.f_install*self.rate_exchange
        else:
            cost_refrigeration = n_ref * \
                (14000*(28.43*P_ref/(T_dispensing_max+273.15)) **
                0.8579+35000*(m_HX/1000))*self.f_install*self.rate_exchange

        return cost_refrigeration

    def cost_booster(self):
        """"Investment of Booster Compressor"""
        booster_compressor_design = self.station.booster_compressor_design()
        n_booster = booster_compressor_design['n_booster']
        P_booster_motor = booster_compressor_design['P_booster_motor']
        cost_booster = n_booster*P_booster_motor*6000*self.f_install*self.rate_exchange

        return cost_booster

    def cost_compressor(self):
        """"Investment of Compressor"""
        compresor_design = self.station.compressor_design()
        n_compressor = compresor_design['n_compressor']
        P_compressor_motor = compresor_design['P_compressor_motor']
        if self.station.dispensing_option == "Cascade Dispensing" and self.station.vehicle_service_pressure == 700:
            cost_compressor = n_compressor*40035 * \
                (P_compressor_motor**0.6038)*self.f_install*self.rate_exchange
        else:
            cost_compressor = n_compressor*40528 * \
                (P_compressor_motor**0.4603)*self.f_install*self.rate_exchange

        return cost_compressor

    def cost_cascade(self):
        """Investment of Cascade"""
        cascade_system = self.station.cascade_calculation()
        capacity_at_pmax_sum = cascade_system['capacity_at_pmax_sum']
        if self.station.dispensing_option == "Booster Compressor Dispensing" or self.station.vehicle_service_pressure == 350:
            cost_cascade = 1100*capacity_at_pmax_sum*self.f_install*self.rate_exchange
        else:
            cost_cascade = 2000*capacity_at_pmax_sum*self.f_install*self.rate_exchange

        return cost_cascade

    def cost_control(self):
        """cost of Overall Control and Safety Equipment"""
        cost_control = 100000*self.rate_exchange
        return cost_control

    def cost_electrical(self):
        """cost of electical"""
        cost_electrical = 80000*self.rate_exchange
        return cost_electrical 

    def investment(self):
        """ Total Capital Investment (€)"""
        investment = {'cost_refrigeration': self.cost_refrigeration(),
                    'cost_booster': self.cost_booster(),
                    'cost_compressor': self.cost_compressor(),
                    'cost_dispenser': self.cost_dispenser(),
                    'cost_cascade': self.cost_cascade(),
                    'cost_lp_storage': self.cost_lp_storage(),
                    'cost_control': self.cost_control(),
                    'cost_electrical': self.cost_electrical()}
        investment_total = 0
        for v in investment.values():
            investment_total += 1.3*v
        investment['investment_total'] = investment_total

        return investment


def calcu_w_spec(senario):
    station = GaseousH2Station(senario)

    compresor_design = station.compressor_design()
    w_sepc_compressor_motor = compresor_design['w_sepc_compressor_motor']

    booster_compresor_design = station.booster_compressor_design()
    w_sepc_booster_motor = booster_compresor_design[
        'w_sepc_booster_motor'] if station.dispensing_option == "Booster Compressor Dispensing" else 0

    refrigeration_design = station.refrigeration_design()
    w_refrigeration = refrigeration_design['w_refrigeration']

    w_spec_total = w_sepc_compressor_motor + w_sepc_booster_motor+w_refrigeration

    w_spec = {'w_sepc_compressor_motor': w_sepc_compressor_motor,
              'w_sepc_booster_motor': w_sepc_booster_motor,
              'w_refrigeration': w_refrigeration,
              'w_spec_total': w_spec_total}

    return w_spec


def calcu_cost_spec(s):
    investment = Investment(s)
    investment_total = investment.investment()['investment_total']
    station = GaseousH2Station(s)

    q = 0.08  # loan rate %
    n = 10  # time [a]
    a = ((1+q)**n)*q/((1+q)**n-1)  # Annuity factor
    annuity = investment_total*a
    demand_year = station.station_size*station.utilization_rate * \
        365  # Sum of hydrogen refueled per year [kg]
    CAPEX = annuity/demand_year

    # fix OPEX
    fix_OPEX = investment_total*0.1/demand_year
    # Variable OPEX
    w_spec_total = calcu_w_spec(s)['w_spec_total']
    price_electric = 0.1509  # €/KWh
    cost_electric = w_spec_total*price_electric

    e_lost_compressor = station.get_input_compressor()[
        'e_lost_compressor']
    price_h2 = 9.5  # €/kg
    cost_lost_compressor = price_h2*e_lost_compressor
    var_OPEX = cost_electric + cost_lost_compressor
    # total cost
    cost_spec_total = CAPEX + fix_OPEX + var_OPEX
    cost_spec = {'CAPEX': CAPEX,
                 'fix_OPEX': fix_OPEX,
                 'var_OPEX': var_OPEX,
                 'cost_spec_total': cost_spec_total}

    return cost_spec
