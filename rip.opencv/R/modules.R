
rip.cv <- new.env(parent = emptyenv())

.initModules <- function()
{
    for (m in c("IO", "core", "transforms", "filter", "feature",
                "photo", "imgproc"))
    {
        rip.cv[[m]] <- Module(m, PACKAGE = "rip.opencv")
    }
    enums <- Module("enums", PACKAGE = "rip.opencv")
    rip.cv[["enums"]] <-
        list(BorderTypes = enums$BorderTypes(),
             DftFlags = enums$DftFlags(),
             InterpolationFlags = enums$InterpolationFlags(),
             ImreadModes = enums$ImreadModes(),
             ColorConversionCodes = unlist(enums$ColorConversionCodes()), # see enums.cpp
             Misc = enums$Misc())
}

## rip.cv$IO <- Module("IO", PACKAGE = "rip.opencv")

