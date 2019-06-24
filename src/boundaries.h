// -----------------------------------------------------------------------------
//
// Copyright (C) The BioDynaMo Project.
// All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
//
// See the LICENSE file distributed with this work for details.
// See the NOTICE file distributed with this work for additional information
// regarding copyright ownership.
//
// -----------------------------------------------------------------------------
#ifndef BOUNDARIES_H_
#define BOUNDARIES_H_

#include "biodynamo.h"

// set numer of simulation steps
const int simulation_steps = 600; // Time between two simulation steps equals: 0.01hours (default)

const double x_range = 150, y_range = 150, z_range = 4500; // set dims of simulation space
//const double simulation_cube_dim = std::max(x_range, y_range);
const double simulation_cube_dim = 4500;  // cube will be 4500x4500x4500

const double z_pos_init = 0 - (simulation_cube_dim/2); // (x/y range + z position) i.e -25
// set boundaries for x/y valuse in the cube
const double x_min = -x_range/2;
const double x_max = x_range/2;
const double y_min = -x_range/2;
const double y_max = x_range/2;


// number of precursor cells
const size_t num_cells = 10;  // number of cells (S1)

namespace bdm {

// 0. Define my custom cell MyCell, which extends Cell by adding extra data
// members: cell_color and can_divide
BDM_SIM_OBJECT(MyCell, Cell) {  // our object extends the Cell object
// create the header with our new data member
BDM_SIM_OBJECT_HEADER(MyCell, Cell, 1, can_divide_, cell_color_);

public:
    MyCellExt() {}
    explicit MyCellExt(const std::array<double, 3>& position) : Base(position) {
        cell_color_[kIdx] = 0;
        can_divide_[kIdx] = true;
    }

   /// If MyCell divides, daughter 2 copies the data members from the mother
    template <typename TMother>
    MyCellExt(const CellDivisionEvent& event, TMother* mother) : Base(event, mother) {
        can_divide_[kIdx] = mother->can_divide_[mother->kIdx];
        cell_color_[kIdx] = mother->cell_color_[mother->kIdx];
    }

    /// If a cell divides, daughter keeps the same state from its mother.
    template <typename TDaughter>
    void EventHandler(const CellDivisionEvent& event, TDaughter* daughter) {
        Base::EventHandler(event, daughter);
    }

    // getter and setter for our new data member
    void SetCanDivide(bool d) { can_divide_[kIdx] = d; }
    bool GetCanDivide() const { return can_divide_[kIdx]; }

    void SetCellColor(int cell_color) { cell_color_[kIdx] = cell_color; }
    int GetCellColor() const { return cell_color_[kIdx]; }

private:
    // declare new data member and define their type
    // private data can only be accessed by public function and not directly
    vec<bool> can_divide_;
    vec<int> cell_color_;
};

// 1. Define growth behaviour
struct GrowthModule : public BaseBiologyModule {
    GrowthModule() : BaseBiologyModule(gAllEventIds) {}

    /// Empty default event constructor, because GrowthModule does not have state.
    template <typename TEvent, typename TBm>
    GrowthModule(const TEvent& event, TBm* other, uint64_t new_oid = 0)
            : BaseBiologyModule(event, other, new_oid) {}

    /// event handler not needed, because Chemotaxis does not have state.
    template <typename T, typename TSimulation = Simulation<>>
    void Run(T* cell) {
        if (cell->GetDiameter() < 8) {
            cell->ChangeVolume(300);
        } else {  //
            if (cell->GetCanDivide()) {
                auto&& daughter = cell->Divide();

                // daughter take the cell_color_ value of her mother
                daughter->SetCellColor(cell->GetCellColor());
                daughter->SetCanDivide(true);  // the daughter will be able to divide
            }
        }

        // bound cells to the x/y range in the cube
        std::array<double, 3> coord = cell->GetPosition();
        bool update_position = false;
        // check for x
        if(coord[0] > x_max - 0.01) {  // if too far in x
            coord[0] = x_max - 0.01;    // reset position to valid value
            update_position = true;
        } else if(coord[0] < x_min + 0.01) {  // if too far in -x
            coord[0] = x_max + 0.01;
            update_position = true;
        }
        // check for y
        if(coord[1] > y_max - 0.01) {  // if too far in y
            coord[1] = y_max - 0.01;    // reset position to valid value
            update_position = true;
        } else if(coord[1] < y_min + 0.01) {  // if too far in -y
            coord[1] = x_max + 0.01;
            update_position = true;
        }
        if (update_position) {
            cell->SetPosition(coord);
        }
    }

    private:
        BDM_CLASS_DEF_NV(GrowthModule, 1);
};

// Define compile time parameter
BDM_CTPARAM() { BDM_CTPARAM_HEADER();
    using SimObjectTypes = CTList<MyCell>;  // use MyCell object
    // Override default BiologyModules for Cell
    BDM_CTPARAM_FOR(bdm, MyCell) { using BiologyModules = CTList<GrowthModule>; };
};

// define actual simulation
inline int Simulate(int argc, const char** argv) {
    // set space parameters of the simulation
    auto set_param = [](auto* param) {
        param->bound_space_ = true;
        param->min_bound_ = -(simulation_cube_dim/2);
        param->max_bound_ = (simulation_cube_dim/2);
        param->run_mechanical_interactions_ = true;
    };

    Simulation<> simulation(argc, argv, set_param);
    auto* rm = simulation.GetResourceManager();  // get pointer to resource manager
    auto* random = simulation.GetRandom();  // get thread of local random number generator.

    double x_coord, y_coord, z_coord;

    // allocate the correct number of cell in our cells structure before
    // cell creation
    rm->template Reserve<MyCell>(num_cells);

    // create 2d Layer of cells
    for (size_t i = 0; i < num_cells; ++i) {
        // create coordinates for cells in 2D plate
        x_coord = random->Uniform(x_min, x_max);
        y_coord = random->Uniform(y_min, y_max);
        z_coord = z_pos_init;

        // creating the cell at position x, y, z
        MyCell cell({x_coord, y_coord, z_coord});
        // set cell parameters
        cell.SetDiameter(6);
        cell.SetAdherence(0.0001);
        cell.SetMass(0.1);
        cell.AddBiologyModule(GrowthModule());
        rm->push_back(cell);// put the created cell in our cells structure
    }


    // 4. Run simulation for N timesteps
    simulation.GetScheduler()->Simulate(simulation_steps);

    std::cout << "Simulation completed successfully!" << std::endl;
    return 0;
}

}  // namespace bdm

#endif  // BOUNDARIES_H_
