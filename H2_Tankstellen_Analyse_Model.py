# -*- coding: utf-8 -*-
"""
Created on Tue Apr 10 19:40:56 2018

@author: Yusheng CHEN
"""
import CoolProp.CoolProp as CP
import numpy as np
import pandas as pd
import math
import json
import GaseoueH2Station as g
import LiquidH2Station as l

def to_excel(obj, writer, sheet_name):
    if isinstance(obj,dict):
        obj = pd.Series(obj,name='python model')
    else:
        pass
    obj.to_excel(writer, sheet_name)


def main():

    # import INPUT DATA
    with open("senario.json", 'r') as s:
        s = json.load(s)

    writer = pd.ExcelWriter('output_10.05.xlsx')
    to_excel(s, writer, 'senario')
    if s['station type'] == 'Gaseous H2 station':
        print('The station ist a gaseous h2 station!')
        station = g.GaseousH2Station(s)
        to_excel(g.calcu_w_spec(s), writer, 'w spec')
        to_excel(g.Investment(s).investment(), writer, 'investment')
        to_excel(g.calcu_cost_spec(s), writer, 'cost spec')        
        to_excel(station.low_pressure_storage_design(), writer, 'low pressure storage')
        to_excel(station.get_input_booster_compressor(), writer, 'input booster compressor')
        to_excel(station.booster_compressor_design(),writer, 'booster compressor')
    else:
        print('The station ist a liquid h2 station!')
        station = l.LiquidH2Station(s)
        to_excel(l.calcu_w_spec(s), writer, 'w spec')
        to_excel(l.Investment(s).investment(), writer, 'investment')
        to_excel(l.calcu_cost_spec(s), writer, 'cost spec')  
        to_excel(station.liquid_cryogenic_tank_design(), writer, 'liquid_cryogenic_tank')
        to_excel(station.pump_design(), writer, ' pump')

    to_excel(station.get_p_vessels(),writer, 'p vessels')
    to_excel(station.get_input_dispenser(), writer, 'input dispenser')
    to_excel(station.refrigeration_design(), writer, 'refirgeration design')
    to_excel(station.get_input_compressor(), writer, 'input compressor')
    to_excel(station.compressor_design(), writer, 'compressor_design')
    to_excel(station.cascade_calculation()['cascade_design'], writer, 'cascade design')
    to_excel(station.cascade_calculation(), writer, 'cascade system')

    writer.save()


if __name__ == "__main__":
    main()
