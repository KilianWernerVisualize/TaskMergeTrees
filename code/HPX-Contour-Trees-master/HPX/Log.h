#pragma once

#include <chrono>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include <unistd.h>

#include <glm/glm.hpp>
#include <hpx/hpx.hpp>

/**
 * @brief Logging utility class.
 */
class Log {
public:
    enum class Level {
        DEBUG_LEVEL,
        INFO_LEVEL,
        WARNING_LEVEL,
        ERROR_LEVEL,
        FILE_LEVEL
    };

public:
    static std::ofstream outfile;
    static hpx::lcos::local::mutex outlock;
    //static hpx::util::high_resolution_timer timer;

    static void init(std::string path){
               //outfile.open(path+std::to_string(hpx::get_locality_id()), std::ios_base::out);
                //timer.restart();
               Log() << "Created log file: " << path+std::to_string(hpx::get_locality_id());
    }

    static void closeFile(){
        //outfile.close();
    }

public:
    Log(Level level = Level::INFO_LEVEL)
        : level(level)
    {
        this->time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    }

    virtual ~Log()
    {
//        if (this->level == Level::DEBUG_LEVEL){
//            if (timer.elapsed() < 600)
//                return;
//        }

        std::ostream* out = &std::cout;
        std::lock_guard<hpx::lcos::local::mutex> lock(Log::outlock);
        if (this->level == Level::FILE_LEVEL){
            out = &Log::outfile;
        }

        std::string s = this->buffer.str();

        if (this->showTime) {
            struct std::tm* ptm = std::localtime(&this->time);
            char buf[80];
            strftime(buf, sizeof(buf), "%T", ptm);
            *out << "[" << std::string(buf) << "] ";
        }

        if (this->showTag.size() > 0) {
            *out << "[" << this->showTag << "] ";
        }

        if (this->level == Level::WARNING_LEVEL)
            *out << "Warning: ";
        else if (this->level == Level::ERROR_LEVEL)
           *out << "Error: ";

        *out << s;

        if (this->newLine)
            *out << std::endl;
        else
            out->flush();
    }

    Log& hideTime()
    {
        this->showTime = false;
        return *this;
    }

    Log& tag(const std::string& tag)
    {
        this->showTag = tag;
        return *this;
    }

    Log& noNewLine()
    {
        this->newLine = false;
        return *this;
    }

    // 2D GLM vectors
    template <typename T, glm::precision P = glm::defaultp>
    Log& operator<<(const glm::tvec2<T, P>& v)
    {
        this->buffer << "(" << v.x << ", " << v.y << ")";
        return *this;
    }

    // 3D GLM vectors
    template <typename T, glm::precision P = glm::defaultp>
    Log& operator<<(const glm::tvec3<T, P>& v)
    {
        this->buffer << "(" << v.x << ", " << v.y << ", " << v.z << ")";
        return *this;
    }

    // 4D GLM vectors
    template <typename T, glm::precision P = glm::defaultp>
    Log& operator<<(const glm::tvec4<T, P>& v)
    {
        this->buffer << "(" << v.x << ", " << v.y << ", " << v.z << ", " << v.w << ")";
        return *this;
    }

    // std::vector
    template <typename T>
    Log& operator<<(const std::vector<T>& v)
    {
        for (int i = 0; i < v.size(); ++i) {
            *this << v[i];
            if (i < v.size() - 1)
                *this << " ";
        }

        return *this;
    }

    // Everything else...
    template <typename T>
    Log& operator<<(const T& t)
    {
        this->buffer << t;
        return *this;
    }

    /**
     * @brief Prints a progress bar
     * @param progress
     */
    void printProgress(float progress)
    {
        int barWidth = 60;

        *this << "[";
        int pos = barWidth * progress;
        for (int i = 0; i < barWidth; ++i) {
            if (i < pos)
                *this << "=";
            else if (i == pos)
                *this << ">";
            else
                *this << " ";
        }
        *this << "] " << int(progress * 100.0) << " %\r";
    }

    void clearLine(int chars = 80)
    {
        this->newLine = false;

        *this << "\r";

        for (int i = 0; i < chars; ++i)
            *this << " ";

        *this << "\r";
    }

private:
    bool newLine = true;

    std::time_t time;
    bool showTime = true;

    std::string showTag;

    Level level;
    std::stringstream buffer;
};

// Error logging
class LogError : public Log {
public:
    LogError()
        : Log(Level::ERROR_LEVEL)
    {
    }
};

// Warning logging
class LogWarning : public Log {
public:
    LogWarning()
        : Log(Level::WARNING_LEVEL)
    {
    }
};

// Info logging
class LogInfo : public Log {
public:
    LogInfo()
        : Log(Level::INFO_LEVEL)
    {
    }
};

class LogFile : public Log {
public:
    LogFile()
        : Log(Level::FILE_LEVEL)
    {
    }
};

// Debug logging
#ifdef ENABLE_DEBUG_LOGGING
class LogDebug : public Log {
public:
    LogDebug()
        : Log(Level::DEBUG_LEVEL)
    {
    }
};
#else
class LogDebug {
public:
    template <typename T>
    LogDebug& operator<<(const T&) { return *this; }

    LogDebug& tag(const std::string&) { return *this; }
};
#endif

// Converts number of bytes to a readable string representation
inline std::string byteString(size_t n, bool addSuffix = true)
{
    std::stringstream s;

    double bytes = (double) n;

    std::string suffix[] = { "B", "KB", "MB", "GB", "TB", "PB", "EB" };

    for (int i = 0; i < 7; ++i) {
        if (bytes < 1024.0) {
            s << bytes;

            if (addSuffix)
                s << " " << suffix[i];

            return s.str();
        }

        bytes /= 1024.0;
    }

    return std::string("wow such many bytes");
}

// Returns the hostname of this machine
inline std::string hostname()
{
    char h[HOST_NAME_MAX];
    gethostname(h, HOST_NAME_MAX);

    return std::string(h);
}
