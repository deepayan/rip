
rip.cv <- new.env(parent = emptyenv())

.initModules <- function()
{
    for (m in c("IO", "core", "transforms", "filter", "feature",
                "photo", "imgproc")) # "xphoto" if available
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
             AdaptiveThresholdTypes = enums$AdaptiveThresholdTypes(),
             DistanceTransformLabelTypes = enums$DistanceTransformLabelTypes(),
             DistanceTransformMasks = enums$DistanceTransformMasks(),
             DistanceTypes = enums$DistanceTypes(),
             FloodFillFlags = enums$FloodFillFlags(),
             GrabCutClasses = enums$GrabCutClasses(),
             GrabCutModes = enums$GrabCutModes(),
             ThresholdTypes = enums$ThresholdTypes(),
             Misc = enums$Misc())
}

## rip.cv$IO <- Module("IO", PACKAGE = "rip.opencv")

