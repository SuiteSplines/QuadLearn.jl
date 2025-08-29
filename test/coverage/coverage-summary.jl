# Coverage summary, printed as "(percentage) covered".
#
# Useful for CI environments that just want a summary (eg a Gitlab setup).

using Coverage

cd(joinpath(@__DIR__, "..", "..")) do
    processed = process_folder()
    covered_lines, total_lines = get_summary(processed)
    percentage = covered_lines / total_lines * 100
    println("($(percentage)%) covered")
    #Codecov.submit_local(processed) ## broken as of 2021-02-23
    # fix: in .gitlab-ci.yml use bash to push coverage.lcov to codecov
    LCOV.writefile("test/coverage/coverage.lcov", processed)
end