module FerSatGui
    using Makie;
    using AbstractPlotting;
    using Observables;
    using DataStructures;

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
    function createCheckBoxGroup(txts::AbstractArray{String, N}, maxOn=UInt64(1)::Unsigned) where {N}
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
            cF = fIndex == -1 ? nothing : functionNamePair[fIndex].func;
            return cF == nothing ? (1:2) : cF(functionNamePair[fIndex].arguments...);
        end
        function safeFunctionCallName(active, index)
            fIndex = length(active) < index ? -1 : active[index];
            cF = fIndex == -1 ? nothing : functionNamePair[fIndex].name;
            return cF == nothing ? "nothing" : cF;
        end
        plots = [lines(lift(safeFunctionCall, checkBoxes.activeBoxIndices, i))[end].parent for i=1:maxPlots];
        texts = [text(lift(safeFunctionCallName, checkBoxes.activeBoxIndices, i), camera=campixel!, raw=true, position=(50, 15))[end].parent for i=1:maxPlots];
        AbstractPlotting.layout([lines(1:10), text("yooo", camera=campixel!, raw=true, position=(50, 15))], 2; sizes=[0.95, 0.05])
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
    function createMultiPlot(funcs::Array{Function, N}, maxPlots) where {N}
        return createMultiPlot([NamedCuryFunction(i, [], " ") for i=funcs], maxPlots);
    end
    function join(multiplot::MultiPlot)
            vbox(hbox([i.scene for i=multiplot.checkBoxes.checkBoxes]...), multiplot.scenes...)
    end
end  # module FerSatGui


#=
createCheckBox("Hey man", false)[1]
boxes = createCheckBoxGroup(["hey" "yo" "what's up" "what's up" "what's up" "what's up"], 4)
on(boxes[2]) do ev
    println(ev)
end
text = AbstractPlotting.text(lift(t -> "curActive: $(string(t))", boxes[2]));
boxes[1] |> i -> hbox([j[1] for j=i]..., text)
=#
