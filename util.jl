# some utilities
function get_arg(arg::String)::Bool
    for a in ARGS
        if a == arg return true end
    end

    return false
end

function get_debug()::Bool
   return get_arg("debug")
end

function get_only_graph()::Bool
   return get_arg("graph")
end

function get_only_data()::Bool
    return get_arg("data")
end
