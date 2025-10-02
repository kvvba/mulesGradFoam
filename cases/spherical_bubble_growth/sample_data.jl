using CSV

tend = 0.6e-3

function extract_lines(input_file, time)
    timesteps = 1
    open("sampled_" * input_file, "w") do output
        open(input_file, "r") do input
            for line in eachline(input)
                if timesteps <= samples && parse(Float64, split.(line)[1]) >= time[timesteps]
                    write(output, line, "\n")
                    timesteps += 1
                end
            end
        end
    end
end


input_list = [
    "05micron.csv",
    "1micron.csv",
    "2micron.csv"
]

print("Number of samples:\n=> ")
samples = parse(Int64, readline())
time_list = [((tend/(samples-1)) * i) for i in 0:(samples-1)]

for i in 1:length(input_list)
    println("Sampling from $(input_list[i]) ...")
    extract_lines(input_list[i], time_list)
end
