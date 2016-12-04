'''
A short script to plot the ELLIE fibre positions on 
a PSUP map.

Author: Ed Leming
Date  : 31/08/2016
'''

import rat
import ROOT

import optparse
import csv
import numpy as np
import sys 
import json
import datetime
import pytz

def read_install_table(fname):
    '''Read in relavent fields from Sofia's install table
    '''
    fibres = []
    nodes = {}
    pmt_hex = {}
    neighbour_hex = {}
    with open(fname, 'rb') as f:
        reader = csv.reader(f)
        for row in reader:
            # Read in the node info - there's one wild one for FT001
            # which isn't a number, just continue past this one.
            try:
                int(row[0])
            except:
                continue
            
            node = int(row[0])
            # Read in fibres, make sure to get the ones wrongly installed!
            if row[3] == '':
                fibre = row[2]
            else:
                fibre = row[3]
            if fibre[:2] == "FT":
                fibres.append("%sA" % fibre)
                fibres.append("%sB" % fibre)
                nodes["%sA" % fibre] = node
                nodes["%sB" % fibre] = node
                # Read in plate positions, make sure to get the ones wrongly installed!
                if row[12] == '':
                    pmt_hex["%sA" % fibre] = row[8]
                    pmt_hex["%sB" % fibre] = row[8]
                    neighbour_hex["%sA" % fibre] = row[9]
                    neighbour_hex["%sB" % fibre] = row[9]
                else: 
                    pmt_hex["%sA" % fibre] = row[12]
                    pmt_hex["%sB" % fibre] = row[12]
                    neighbour_hex["%sA" % fibre] = row[13]
                    neighbour_hex["%sB" % fibre] = row[13]
            else: 
                fibres.append(fibre)
                nodes[fibre] = node
                # Read in plate positions, make sure to get the ones wrongly installed!
                if row[12] == '':
                    pmt_hex[fibre] = row[8]
                    neighbour_hex[fibre] = row[9]
                else: 
                    pmt_hex[fibre] = row[12]
                    neighbour_hex[fibre] = row[13]

    return nodes, fibres, pmt_hex, neighbour_hex

def read_PMT_coordinates(fname):
    '''Get global positions of PMT Hex cells
    '''
    cells = {}
    with open(fname, 'rb') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if row == []:
                continue
            row = filter(None, row)
            try:
                cells[row[1]] = [float(row[3]), float(row[4]), float(row[5])]
            except:
                raise
    return cells
    
def get_pmt_coordinates(host_hex, neighbour_hex, fname):
    '''Use the PMT hex cell file to find PMT positions and return
    the relavent three vectors for the plate position calculation
    '''
    # Some brute force for correcting fields read from .csv
    if "none" in neighbour_hex:
        neighbour_hex = host_hex

    panel_pmts = {}
    with open(fname, 'rb') as f:
        reader = csv.reader(f, delimiter='\t')
        next(reader, None)  # Skip header line
        for row in reader:
           # Skip over empty lines
            if row == []:
                continue
            # Remove empty strings
            row = filter(None, row)
            # Get three vectors for host and neighbour cells + record the hex_panel name
            if row[1][:4] == host_hex:
                host_vec = ROOT.TVector3(float(row[3]), float(row[4]), float(row[5]))
                panel_name = row[6][:4]
            if row[1][:4] == neighbour_hex:
                neighbour_vec = ROOT.TVector3(float(row[3]), float(row[4]), float(row[5]))
        # If no neighbour listed then panel is directly above host
        try:
            neighbour_vec
        except:
            neighbour_vec = host_vec

        # Now find all cells in that panel!
        f.seek(0)
        next(reader, None) # Skip header line
        for row in reader:
            # Skip over empty lines
            if row == []:
                continue
            # Remove empty strings
            row = filter(None, row)
            if row[6][:4] == panel_name:
                panel_pmts[int(row[6][-2:])] = [float(row[3]), float(row[4]), float(row[5])]

    # Find central PMT in this hex panel.
    if len(panel_pmts) > 9:
        central_pmt = panel_pmts[10]
    else:
        central_pmt = panel_pmts[4]
    central_vec = ROOT.TVector3(central_pmt[0], central_pmt[1], central_pmt[2])

    # Extract a panel number for querying the .ratdb
    panel_number =  int(panel_name[1:])
    
    # return
    return central_vec, host_vec, neighbour_vec, panel_number

def calc_fibre_placement_TELLIE(host_hex, neighbour_hex, fname):
    '''Use the install hex and neighbour hex to calculate the plate install position
    '''
    # Define hex cell & plate parameters
    hex_size = 17.05                             # [cm]
    plate_height = 0.3 + 0.22                    # [cm]
    hex_radius = float(hex_size*np.cos(30*(np.pi/180))) # [cm]
    fibre_separation = 0.51                      # [cm]
    
    vec_from_centre, host_vec, neighbour_vec, panel_number = get_pmt_coordinates(host_hex, neighbour_hex, fname)

    ##########################################
    # Points to the centre of the detctor!
    ##########################################
    # You can use the PMT at the centre of the panel
    # to give a pointing vector,
    # vec_to_centre = -vec_from_centre
    ##########################################
    # Alternatively, use pointing vectors in PANELINFO.ratdb. 
    panel_info = ROOT.RAT.DB.Get().GetLink("PANELINFO", "")
    panels = panel_info.GetIArray("panel_number")
    panel_list = []
    for i, p in enumerate(panels):
        panel_list.append(panels[i])
    db_index = panel_list.index(panel_number)
    vec_to_centre = ROOT.TVector3(panel_info.GetDArray("u")[db_index],
                                  panel_info.GetDArray("v")[db_index],
                                  panel_info.GetDArray("w")[db_index])
    ##########################################
    vec_to_centre.SetMag(1.) 
    
    # Unit vector pointing from central hex to neighbour hex
    cent_to_neig_unit = neighbour_vec - vec_from_centre
    cent_to_neig_unit.SetMag(1.0)

    # Vector to mounted plate
    plate_vec = host_vec + cent_to_neig_unit*(hex_radius + plate_height)
    
    # Correct for 3D curve of PSUP
    curve_correction = (plate_vec - host_vec).Cross(vec_to_centre)
    curve_correction.SetMag(1.)

    # A & B positions
    positionA = plate_vec + curve_correction*fibre_separation
    positionB = plate_vec - curve_correction*fibre_separation

    ################################################################
    # Some checks lifted stright from Sofia's code (docdb 1730)
    #
    ################################################################
    #da=float(np.sin(positionA.Angle(vec_to_centre))*positionA.Mag())
    #db=float(np.sin(positionB.Angle(vec_to_centre))*positionB.Mag())
    
    #aa=cent_to_neig_unit.Angle(vec_to_centre)*(180./np.pi)
    #ab=(positionA-positionB).Angle(vec_to_centre)*(180./np.pi)

    #if(da>500. or db>500. or np.abs(aa-90.)>0.5 or np.abs(ab-90.)>0.5):
    #    print "dist=", (positionA-positionB).Mag(), " in R=", positionA.Mag()-positionB.Mag()
    #    print "to center=", da, " and ",db
    #    print "Angle cor x dir=", aa, " and poA-poB x dir="

    return positionA, positionB

def calc_fibre_placement_SMELLIE(host_hex, neighbour_hex, fname):
    '''Use the install hex and neighbour hex to calculate the plate install position
    '''
    # Define hex cell & plate parameters
    hex_size = 17.05                             # [cm]
    plate_height = 0.15                          # [cm]
    hex_radius = float(hex_size*np.cos(30*(np.pi/180))) # [cm]
    abrplace1 = 0.1586+0.05+0.0779               # [cm]
    abrplace2 = 0.1237+0.05+0.0735               # [cm]
    # fibres come past the end of the plate
    protrude = 0.025                             # [cm]
    protrude10 = 0.025 - 0.014                     # [cm]
    protrude20 = 0.025 - 0.027                     # [cm]    

    vec_from_centre, host_vec, neighbour_vec, panel_number = get_pmt_coordinates(host_hex, neighbour_hex, fname)

    ##########################################
    # Points to the centre of the detctor!
    ##########################################
    # You can use the PMT at the centre of the panel
    # to give a pointing vector,
    #vec_to_centre = -vec_from_centre
    ##########################################
    # Alternatively, use pointing vectors in PANELINFO.ratdb. 
    panel_info = ROOT.RAT.DB.Get().GetLink("PANELINFO", "")
    panels = panel_info.GetIArray("panel_number")
    panel_list = []
    for i, p in enumerate(panels):
        panel_list.append(panels[i])
    db_index = panel_list.index(panel_number)
    vec_to_centre = ROOT.TVector3(panel_info.GetDArray("u")[db_index],
                                  panel_info.GetDArray("v")[db_index],
                                  panel_info.GetDArray("w")[db_index])
    ##########################################
    vec_to_centre.SetMag(1.) 
    
    # Unit vector pointing from central hex to neighbour hex
    cent_to_neig_unit = neighbour_vec - vec_from_centre
    cent_to_neig_unit.SetMag(1.0)

    # Vector to mounted plate
    plate_vec = host_vec + cent_to_neig_unit*(hex_radius + plate_height)
    
    # Correct for 3D curve of PSUP
    curve_correction = (plate_vec - host_vec).Cross(vec_to_centre)
    curve_correction.SetMag(1.)

    # Protrusion corrections
    protrusion_correction_0 = vec_to_centre*protrude
    protrusion_correction_10 = vec_to_centre*protrude10
    protrusion_correction_20 = vec_to_centre*protrude20

    # 0, 10 and 20 deg fibre positions
    position_0 = plate_vec + protrusion_correction_0
    position_10 = plate_vec + curve_correction*abrplace1 + protrusion_correction_10
    position_20 = plate_vec - curve_correction*abrplace2 + protrusion_correction_20

    # 0, 10 and 20 deg fibre pointing angles
    pointing_0 = ROOT.TVector3(vec_to_centre)
    pointing_10 = ROOT.TVector3(vec_to_centre)
    pointing_20 = ROOT.TVector3(vec_to_centre)
    pointing_10.Rotate(10.*(np.pi/180), cent_to_neig_unit)
    pointing_20.Rotate(-20.*(np.pi/180), cent_to_neig_unit)

    return position_0, position_10, position_20, pointing_0, pointing_10, pointing_20


def get_pointing_angle(fibre):
    '''Get pointing angle of a fibre according to the info in the current
    database entry for that fibre
    '''
    try:
        entry = ROOT.RAT.DB.Get().GetLink("FIBRE", fibre)
        # Create vectors pointing to origin at centre of detector
        central_vector = ROOT.TVector3(entry.GetD("x"), entry.GetD("y"), entry.GetD("z"))
        pointing_vector = ROOT.TVector3(entry.GetD("u"), entry.GetD("v"), entry.GetD("w"))
    except Exception as e:
        print "Refernce to fibre %s does not exist in local ratdb" % fibre
        #raise e
    # Return the angle between the two vectors
    return 180 - central_vector.Angle(pointing_vector)*(180/np.pi)
        
def compare_position_calculations_TELLIE(fibres, pmt_hex, neighbour_hex, fname):
    '''Compare the x,y,z positions from our calculations with those in 
    the database
    '''
    count = 0
    for fibre in fibres:
        if fibre[:2] == 'FT':
            table = ROOT.RAT.DB.Get().GetDefaultTable("FIBRE", fibre)
            posA, posB = calc_fibre_placement_TELLIE(pmt_hex[fibre], neighbour_hex[fibre], fname)
            try:
                table.GetD("x")
            except:
                print "###############################"
                print "Coundn't access x, y, z variables in database for fibre %s, skipping" % fibre
                continue        
            
            # Check for differences
            if fibre[-1] == "B":
                diff = table.GetD("x")/10. - posB.x()
            else:
                diff = table.GetD("x")/10. - posA.x()
                
            if np.abs(diff) > 1:
                print "###############################"
                print "%s:" % fibre
                print "Host: %s\nNeighbour: %s" % (pmt_hex[fibre], neighbour_hex[fibre])
                print "Database: \t%1.3f, %1.3f, %1.3f" % (table.GetD("x")/10., table.GetD("y")/10., table.GetD("z")/10.)
                if fibre[-1] == "B":
                    print "Calculated :\t%1.3f, %1.3f, %1.3f" % (posB.x(), posB.y(), posB.z())
                else:
                    print "Calculated :\t%1.3f, %1.3f, %1.3f" % (posA.x(), posA.y(), posA.z())
                        
def compare_position_calculations_SMELLIE(fibres, pmt_hex, neighbour_hex, fname):
    '''Compare the x,y,z positions from our calculations with those in 
    the database
    '''
    count = 0
    for fibre in fibres:
        if fibre[:2] != 'FT':
            print fibre
            table = ROOT.RAT.DB.Get().GetDefaultTable("FIBRE", fibre)
            pos_0, pos_10, pos_20, point_0, point_10, point_20 = calc_fibre_placement_SMELLIE(pmt_hex[fibre], neighbour_hex[fibre], fname)
            try:
                table.GetD("x")
            except:
                print "###############################"
                print "Coundn't access x, y, z variables in database for fibre %s, skipping" % fibre
                continue        
        
            # Check for differences
            database_point = ROOT.TVector3(table.GetD("u"), table.GetD("v"), table.GetD("w"))
            if fibre[2] == "0":
                pos_diff = table.GetD("x")/10. - pos_0.x()
                calc_point = ROOT.TVector3(point_0.x(), point_0.y(), point_0.z())
            elif fibre[2] == "1":
                pos_diff = table.GetD("x")/10. - pos_10.x()
                calc_point = ROOT.TVector3(point_10.x(), point_10.y(), point_10.z())
            elif fibre[2] == "2":
                pos_diff = table.GetD("x")/10. - pos_20.x()
                calc_point = ROOT.TVector3(point_20.x(), point_20.y(), point_20.z())
            
            point_diff = database_point.Angle(calc_point)*(180/np.pi)

            if np.abs(pos_diff) > 1 or np.abs(point_diff) > 1:
                print "###############################"
                print "%s:" % fibre
                print "Host: %s\nNeighbour: %s" % (pmt_hex[fibre], neighbour_hex[fibre])
                print "Database pos: \t%1.3f, %1.3f, %1.3f" % (table.GetD("x")/10., table.GetD("y")/10., table.GetD("z")/10.)
                if fibre[2] == "0":
                    print "Calculated pos:\t%1.3f, %1.3f, %1.3f" % (pos_0.x(), pos_0.y(), pos_0.z())
                elif fibre[2] == "1":
                    print "Calculated pos:\t%1.3f, %1.3f, %1.3f" % (pos_10.x(), pos_10.y(), pos_10.z())
                if fibre[2] == "2":
                    print "Calculated pos:\t%1.3f, %1.3f, %1.3f" % (pos_20.x(), pos_20.y(), pos_20.z())

                print "Difference in angles : %1.2f degrees" % (point_diff)
                print "Database point: \t%1.3f, %1.3f, %1.3f" % (table.GetD("u"), table.GetD("v"), table.GetD("w"))
                if fibre[2] == "0":
                    print "Calculated point:\t%1.3f, %1.3f, %1.3f" % (point_0.x(), point_0.y(), point_0.z())
                elif fibre[2] == "1":
                    print "Calculated point:\t%1.3f, %1.3f, %1.3f" % (point_10.x(), point_10.y(), point_10.z())
                if fibre[2] == "2":
                    print "Calculated point:\t%1.3f, %1.3f, %1.3f" % (point_20.x(), point_20.y(), point_20.z())


def make_new_db_files(fibres):
    '''Make new database files for jose
    '''
    for fibre in fibres:
        if fibre[:2] == "FT":
            angle = 0
        else:
            try:
                angle = int(round(get_pointing_angle(fibre), -1))
            except:
                "Fibre: %s does not exist in database" % fibre
                continue
        #print fibre, nodes[fibre], angle
        table = ROOT.RAT.DB.Get().GetDefaultTable("FIBRE", fibre)
        table.SetI("psup_panel", nodes[fibre])
        table.SetI("pointing_angle", angle)
        table.SetI("fibre_status", 0)
        table.SaveAs("./new_tables/%s.json" % fibre)

def make_pannel_to_fibre_map(fibres,nodes):
    '''Make new node->fibre mapping file for
    pushing to the telliedb. This file is 
    used by orca to map fibres to pannels.
    '''
    timestamp = datetime.datetime.now(pytz.timezone('US/Eastern')).isoformat()

    new_doc = {}
    new_doc["type"] = "PANEL_FIBRE_MAPPING"
    new_doc["comment"] = ""
    new_doc["index"] = ""
    new_doc["version"] = 0
    new_doc["pass"] = 0
    new_doc["first_valid"] = 10870
    new_doc["run_range"] = [10870, 2147483647]
    new_doc["timestamp"] = timestamp
    node_summary = [dict() for x in range(92)]
    for fibre, node in nodes.iteritems():
        if fibre[:2] == "FT":
            node_summary[node-1][fibre] = 0
    for i, node in enumerate(node_summary):
        new_doc["panel_%i" % (i + 1)] = node 
        print "panel_%i" % (i + 1), new_doc["panel_%i" % (i + 1)]
    # Save to file
    with open('node_fibre_mapping.json', 'w') as fp:
        json.dump(new_doc, fp)


if __name__ == "__main__":
    parser = optparse.OptionParser()
    parser.add_option("-i",dest="install",
                      help="Path to csv file detailing the current install status",
                      default="./ellie-fibres/install_table.csv")
    parser.add_option("-p",dest="pmt",
                      help="Path to csv file detailing the global co-ordinates of the PMTs",
                      default="PostRotation.txt")
    (options,args) = parser.parse_args()

    # Reset all root stuff
    ROOT.gROOT.Reset()

    # In in fibre install locations & hex cell positions
    nodes, fibres, pmt_hex, neighbour_hex = read_install_table(options.install)

    # Load rat defualts for fibre pos stuff
    ROOT.RAT.DB.Get().LoadDefaults()

    # make new database files for jose
    #make_new_db_files(fibres)

    # Compare position calculations
    #compare_position_calculations_TELLIE(fibres, pmt_hex, neighbour_hex, options.pmt)
    #compare_position_calculations_SMELLIE(fibres, pmt_hex, neighbour_hex, options.pmt)

    make_pannel_to_fibre_map(fibres, nodes)
