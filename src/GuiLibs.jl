module FerSatGui
    using Makie;
    using AbstractPlotting;
    using Observables;
    using DataStructures;
    using GeometryTypes;

    struct CheckBox
        scene::Scene
        checked::Observable{Bool}
        button::Button
    end
    struct CheckBoxGroup{N, M}
        checkBoxes::AbstractArray{CheckBox, N}
        activeBoxIndices::Observable{Array{Unsigned, M}}
    end
    struct NamedCuryFunction{T, N}
        func::Any
        arguments::Array{T, N}
        name::String
    end
    struct MultiPlot
        scenes::AbstractArray{Scene, 1}
        checkBoxes::CheckBoxGroup
        maxActive::Unsigned
    end
    struct NamedObsevableArray{T, N} 
        array::Observable{<:AbstractArray{T, N}}
        name::String
    end
    function Base.:*(rec::GeometryTypes.HyperRectangle{N, T}, fac) where {N, T}
        GeometryTypes.HyperRectangle(rec.origin..., (rec.widths .* fac)...)
    end
    function Base.:+(rec::GeometryTypes.HyperRectangle{N, T}, pls) where {N, T}
        GeometryTypes.HyperRectangle((rec.origin + pls)..., rec.widths...)
    end
    function createNormalText(_text, pos)
        text(_text; camera=campixel!, raw=true, position=pos)
    end
    function createCheckBox(txt="checkbox"::String, checked=false::Bool)
        s = Scene(camera=AbstractPlotting.campixel!);
        checkedNode = Node(checked);
        textLabel = lift(c -> "$(txt) $(c ? "|✅|" : "|❌|")", checkedNode);
        b = button!(s, textLabel, raw=true, align=(:bottom, :right), position=(10, 10))[end];
        on(b[:clicks]) do ev
            if !isempty(ev)
                push!(checkedNode, !to_value(checkedNode));
            end
        end
        CheckBox(s, checkedNode, b)
    end

    function doItRight(myQ, checkedBoxes, curi, maxOn, curActive)

        on(checkedBoxes[curi].checked) do ev
            if ev
                if !any(myQ .|> j -> j == curi)
                    enqueue!(myQ, curi);
                    push!(curActive, [myQ...]);
                end
                if maxOn < length(myQ)
                    id = dequeue!(myQ);
                    push!(checkedBoxes[id].checked, false);
                    push!(curActive, [myQ...]);
                end
            elseif any(myQ .|> j -> j == curi)
                tempQ = Queue{Integer}();
                while length(myQ)>0
                    if front(myQ) != curi
                        enqueue!(tempQ, dequeue!(myQ));
                    else
                        dequeue!(myQ);
                    end
                end
                while length(tempQ) > 0
                    enqueue!(myQ, dequeue!(tempQ));
                end
                push!(curActive, [myQ...]);
            end

        end
    end

    function createCheckBoxGroup(txts::AbstractArray{String, N}, maxOn=UInt64(1)::Unsigned)::CheckBoxGroup where {N}
        checkedBoxes = createCheckBox.(txts);
        myQ = Queue{Integer}();
        curActive = Node{Array{Unsigned, 1}}([])
        for i=1:length(checkedBoxes)
            doItRight(myQ, checkedBoxes, i, maxOn, curActive);
        end
        CheckBoxGroup(checkedBoxes, curActive)
    end

    function createMultiPlot(functionNamePair::Array{NamedCuryFunction{T, N}, M}, maxPlots=UInt64(2)::Unsigned) where {T, N, M}
        checkBoxes = createCheckBoxGroup([i.name for i=functionNamePair], maxPlots);
        function safeFunctionCall(active, index)
            fIndex = length(active) < index ? -1 : active[index];
            cF = fIndex == -1 ? (1:2) : functionNamePair[fIndex].func(functionNamePair[fIndex].arguments...);
            return cF;
        end
        function safeFunctionCallName(active, index)
            fIndex = length(active) < index ? -1 : active[index];
            cF = fIndex == -1 ? "nothing" : functionNamePair[fIndex].name;
            return cF;
        end

        plots = [lines(lift(safeFunctionCall, checkBoxes.activeBoxIndices, i))[end].parent for i=1:maxPlots];
        texts = [text(lift(safeFunctionCallName, checkBoxes.activeBoxIndices, i), camera=campixel!, raw=true, position=(50, 15))[end].parent for i=1:maxPlots];
        # AbstractPlotting.layout([lines(1:10), text("yooo", camera=campixel!, raw=true, position=(50, 15))], 2; sizes=[0.95, 0.05])
        on(checkBoxes.activeBoxIndices) do ev
            for i=1:length(ev)
                if (length(plots) >= i)
                    AbstractPlotting.update_limits!(plots[i]);
                    AbstractPlotting.update_cam!(plots[i]);
                end
            end
        end
        scenes = Array{Scene, 1}();
        for i=1:length(plots)
            push!(scenes, AbstractPlotting.layout([plots[i], texts[i]], 2; sizes=[0.95, 0.05]));
        end
        MultiPlot(scenes, checkBoxes, maxPlots);
    end
    function doIt(lifts, i, liftsToNodes, nodes, myLines, textNodes) 
        namedOA = lifts[i];
        on(namedOA.array) do ev
            if haskey(liftsToNodes, i)
                nodeToUpdate = liftsToNodes[i];
                push!(nodes[nodeToUpdate], ev);
                # AbstractPlotting.update_limits!(myLines[nodeToUpdate]);
                # AbstractPlotting.update_cam!(myLines[nodeToUpdate]);
                # AbstractPlotting.update_cam!(myLines[nodeToUpdate], area) = update_cam!(scene, cameracontrols(scene), area)
                push!(textNodes[nodeToUpdate], lifts[i].name)
            end
        end
    end
    function createMultiPlot(lifts::Array{NamedObsevableArray{T, N}, M}, maxPlots) where {T1 <:Real, T <: Union{T1, Tuple{T1, T1}}, N, M} 
        checkboxes::CheckBoxGroup = createCheckBoxGroup([n.name for n=lifts], maxPlots);
        nodes = T  <: Tuple ? [Node{Array{T, N}}([(0, 0), (0, 10), (0, 0), (10, 0)]) for i=1:maxPlots] : [Node{Array{T, N}}([(T1(1):T1(2))...]) for i=1:maxPlots];
        textNodes = [Node("$i") for i=1:maxPlots];
        liftsToNodes = Dict();
        myLines = [lines(n)[end].parent for n=nodes];
        myTexts = [text(n, camera=campixel!, raw=true, position=(50, 15))[end].parent for n=textNodes];
        for i=1:length(lifts)
            doIt(lifts, i, liftsToNodes, nodes, myLines, textNodes);
        end
        on(checkboxes.activeBoxIndices) do ev
            for i=1:length(ev)
                if (length(nodes) >= i)
                    liftIndex = ev[i];
                    liftArray = to_value(lifts[liftIndex].array);
                    push!(nodes[i], liftArray);
                    AbstractPlotting.update_limits!(myLines[i]);
                    AbstractPlotting.scale_scene!(myLines[i]);
                    AbstractPlotting.update_cam!(myLines[i]);
                    # AbstractPlotting.update_limits!(myLines[i], myLines[i][:limits][], myLines[i][:padding][])
                    

                    limits = AbstractPlotting.limits(myLines[i])[];
                    translate = limits.widths .* 0.2;
                    limits = GeometryTypes.HyperRectangle((limits.origin .- translate)..., limits.widths...) * 1.2;

                    AbstractPlotting.update_cam!(myLines[i], limits);
                    push!(textNodes[i], lifts[liftIndex].name)
                    filter!(j -> j[2] != i, liftsToNodes);
                    setindex!(liftsToNodes, i, ev[i]);
                end
            end
        end
        scenes = Array{Scene, 1}();
        for i=1:length(myLines)
            push!(scenes, AbstractPlotting.layout([myLines[i], myTexts[i]], 2; sizes=[0.95, 0.05]));
        end
        return MultiPlot(scenes, checkboxes, maxPlots)
    end

    function createMultiPlot(funcs::Array{Function, N}, maxPlots) where {N}
        return createMultiPlot([NamedCuryFunction(i, [], " ") for i=funcs], maxPlots);
    end

    function join(multiplot::MultiPlot)
            vbox(hbox([i.scene for i=multiplot.checkBoxes.checkBoxes]...), multiplot.scenes...)
    end
    function join(multiplot::MultiPlot, additionalControll::Array{T,N}) where {T<:AbstractPlotting.Transformable, N}
        sizes = [0.25, [0.75/length(multiplot.scenes) for i=multiplot.scenes]...];
        AbstractPlotting.layout([hbox([i.scene for i=multiplot.checkBoxes.checkBoxes]..., additionalControll...), multiplot.scenes...], 1; sizes=sizes)
    end
    function test()
        nodes = [NamedObsevableArray(Node([(1:i)...]), "Do $i") for i=10:10:100];
        mp = createMultiPlot(nodes, 2);
        return join(mp);
    end
end 
