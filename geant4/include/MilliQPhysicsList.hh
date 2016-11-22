
#ifndef MilliQPhysicsList_h
#define MilliQPhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>


class MilliQPhysicsList: public G4VModularPhysicsList
{
  public:

    MilliQPhysicsList(const boost::property_tree::ptree pt);
    virtual ~MilliQPhysicsList();

  public:

    virtual void SetCuts();

};

#endif
