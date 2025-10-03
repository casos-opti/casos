% SPDX-FileCopyrightText: 2023 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

classdef FunctionArgumentType
% Possible types of arguments for functions.

    properties
        classname;
    end

    enumeration
        % vectors & matrices
        DM ('casadi.DM')
        SX ('casadi.SX')
        MX ('casadi.MX') 

        % polynomials
        PD ('casos.PD')
        PS ('casos.PS')
    end

    methods
        function type = FunctionArgumentType(cls)
            % Create new function argument type.
            type.classname = cls;
        end

        function tf = isOfType(type,var)
            % Check if variable is of given type.
            tf = isa(var, type.classname);
        end
    end
end
